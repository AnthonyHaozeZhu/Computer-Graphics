#include "server/Server.hpp"

#include "PathPhotonMapping.hpp"

#include "VertexTransformer.hpp"
#include "intersections/intersections.hpp"

#include "glm/gtc/matrix_transform.hpp"
#include "Onb.hpp"
#include <omp.h>

#define KDTREE

namespace PathPhotonMapping
{
    RGB PathPhotoMappingRenderer::gamma(const RGB& rgb) {
        return glm::sqrt(rgb);
    }
    int index = 0;

    Vec3 randomDir(const Vec3& normal) {
        Vec3 random = defaultSamplerInstance<HemiSphere>().sample3d();
        Onb onb{ normal };
        Vec3 direction = glm::normalize(onb.local(random));
        return direction;
    }
    void PathPhotoMappingRenderer::pthotoMap() {
        for (int i = 0; i < photonNum; i++) {
            for (int j = 0; j < scene.areaLightBuffer.size(); j++) {
                auto p1 = (rand() % 1000) / 1000.f;
                auto p2 = (rand() % 1000) / 1000.f;
                Vec3 random = randomDir(scene.areaLightBuffer[j].n);
                auto origin = scene.areaLightBuffer[j].position
                    + p1 * scene.areaLightBuffer[j].u
                    + p2 * scene.areaLightBuffer[j].v;
                Ray ray{origin, random};

                Vec3 r = (scene.areaLightBuffer[j].radiance / ((1.f / scene.areaLightBuffer[j].A) * PI));
                pthotoMap(ray, r, 0);
            }
        }
        char log[512];
        sprintf(log, "photon num: %u", photons.size());
        getServer().logger.log(log);
        photon_tree.setPoints(photons);
    }

    void PathPhotoMappingRenderer::pthotoMap(Ray& ray, Vec3 r, int currDepth) {
        auto hitObject = closestHitObject(ray);
        if (!hitObject) return;
        else {
            auto mtlHandle = hitObject->material;
            auto scattered = shaderPrograms[mtlHandle.index()]->shade(ray, hitObject->hitPoint, hitObject->normal);
            auto attenuation = scattered.attenuation;
            auto emitted = scattered.emitted;
            Ray next = scattered.ray;

            if(scene.materials[mtlHandle.index()].hasProperty("diffuseColor")){
                photon p{ hitObject->hitPoint, r, ray, next, hitObject->normal};
                photons.push_back(p);
                if (currDepth >= depth || rand() % 1000 < 200) return;
                float cos_0 = abs(glm::dot(hitObject->normal, ray.direction));
                pthotoMap(next, attenuation * r * cos_0 / (scattered.pdf * 0.8f) , currDepth + 1);
                
            }
            else if (scene.materials[mtlHandle.index()].hasProperty("ior")) {
                if (currDepth >= depth || rand() % 1000 < 200) return;
                pthotoMap(next, attenuation * r / (scattered.pdf * 0.8f) , currDepth + 1);
                next = scattered.refraction_ray;
                pthotoMap(next, scattered.refraction * r / (scattered.pdf * 0.8f) , currDepth + 1);
            }
        }
    }

    void PathPhotoMappingRenderer::renderTask(RGBA* pixels, int width, int height) {
        #pragma omp parallel for
        for(int i=0; i<height; i++) {
            for (int j=0; j<width; j++) {
                Vec3 color{0, 0, 0};
                for (int k=0; k < samples; k++) {
                    auto r = defaultSamplerInstance<UniformInSquare>().sample2d();
                    float rx = r.x;
                    float ry = r.y;
                    float x = (float(j)+rx)/float(width);
                    float y = (float(i)+ry)/float(height);
                    auto ray = camera.shoot(x, y);
                    color += (trace_direct(ray, 0) + trace_indirect(ray, 0)) / 2.f;
                }
                color /= samples;
                color = gamma(color);
                pixels[(height-i-1)*width+j] = {color, 1};
            }
        }
    }

    Vec3 outerProduct(const Vec3& a, const Vec3& b) {
        auto x = a.y * b.z - b.y * a.z;
        auto y = -(a.x * b.z - b.x * a.z);
        auto z = a.x * b.y - b.x * a.y;
        return Vec3{ x,y,z };
    }

    auto PathPhotoMappingRenderer::render() -> RenderResult {
        // shaders
        shaderPrograms.clear();
        ShaderCreator shaderCreator{};
        for (auto& m : scene.materials) {
            shaderPrograms.push_back(shaderCreator.create(m, scene.textures));
        }

        RGBA* pixels = new RGBA[width*height]{};

        // 局部坐标转换成世界坐标
        VertexTransformer vertexTransformer{};
        vertexTransformer.exec(spScene);
        tree.Insert(scene.sphereBuffer, scene.triangleBuffer);
        for (int i = 0; i < scene.areaLightBuffer.size(); i++) {
            auto length_u = glm::length(scene.areaLightBuffer[i].u);
            auto length_v = glm::length(scene.areaLightBuffer[i].v);
            auto norm_u = glm::normalize(scene.areaLightBuffer[i].u);
            auto norm_v = glm::normalize(scene.areaLightBuffer[i].v);
            auto sin_0 = std::pow(1 - std::pow(glm::dot(norm_u, norm_v), 2), 0.5);
            float A = length_u * length_v * sin_0;
            auto n = outerProduct(norm_u, norm_v);
            scene.areaLightBuffer[i].A = A;
            scene.areaLightBuffer[i].n = n;
        }
        pthotoMap();
        PathPhotoMappingRenderer::renderTask(pixels, width, height);
        getServer().logger.log("Done...");
        return {pixels, width, height};
    }

    void PathPhotoMappingRenderer::release(const RenderResult& r) {
        auto [p, w, h] = r;
        delete[] p;
    }

    HitRecord PathPhotoMappingRenderer::closestHitObject(const Ray& r) {
        HitRecord closestHit = nullopt;
        float closest = FLOAT_INF;
#ifndef KDTREE
        for (auto& s : scene.sphereBuffer) {
            auto hitRecord = Intersection::xSphere(r, s, 0.01, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        for (auto& t : scene.triangleBuffer) {
            auto hitRecord = Intersection::xTriangle(r, t, 0.01, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }
#endif // KDTREE
#ifdef KDTREE
        auto result = tree.find_Node(r);

        for (auto& n : result) {
            for (auto& b : n->Box_list) {
                HitRecord hitRecord;
                if (b->type == SPHERE) hitRecord = Intersection::xSphere(r, *(b->sp), 0.01, closest);
                else if (b->type == TRIANGLE) hitRecord = Intersection::xTriangle(r, *(b->tr), 0.01, closest);
                if (hitRecord && hitRecord->t < closest) {
                    closest = hitRecord->t;
                    closestHit = hitRecord;
                }
            }
        }
#endif // KDTREE
        for (auto& p : scene.planeBuffer) {
            auto hitRecord = Intersection::xPlane(r, p, 0.000001, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        return closestHit; 
    }
    
    tuple<float, Vec3> PathPhotoMappingRenderer::closestHitLight(const Ray& r) {
        Vec3 v = {};
        HitRecord closest = getHitRecord(FLOAT_INF, {}, {}, {});
        for (auto& a : scene.areaLightBuffer) {
            auto hitRecord = Intersection::xAreaLight(r, a, 0.000001, closest->t);
            if (hitRecord && closest->t > hitRecord->t) {
                closest = hitRecord;
                v = a.radiance;
            }
        }
        return { closest->t, v };
    }

    RGB PathPhotoMappingRenderer::trace_direct(const Ray& r, int currDepth) {
        if (currDepth > depth && rand() % 100 < 20) return scene.ambient.constant;
        auto hitObject = closestHitObject(r);
        auto [ t, emitted ] = closestHitLight(r);
        // hit object
        if (hitObject && hitObject->t < t) {
            auto mtlHandle = hitObject->material;
            auto scattered = shaderPrograms[mtlHandle.index()]->shade(r, hitObject->hitPoint, hitObject->normal);
            auto scatteredRay = scattered.ray;
            auto attenuation = scattered.attenuation;
            auto emitted = scattered.emitted;

            RGB next = Vec3(0.f);
            if (scene.materials[mtlHandle.index()].hasProperty("ior")) {
                if (currDepth > depth && rand() % 1000 < 200) next = (attenuation + scattered.refraction) * scene.ambient.constant;
                else {
                    auto reflex = trace_direct(scatteredRay, currDepth + 1);
                    RGB refraction = Vec3(0.f);
                    if (scattered.refraction_ray.direction != Vec3(0.f))  refraction = trace_direct(scattered.refraction_ray, currDepth + 1);
                    next = attenuation * reflex + refraction * scattered.refraction;
                }
            }
            float pdf = scattered.pdf;

            /**
             * emitted      - Le(p, w_0)
             * next         - Li(p, w_i)
             * n_dot_in     - cos<n, w_i>
             * atteunation  - BRDF
             * pdf          - p(w)
            **/
            auto Lo = (next / pdf)  / 0.8f;
            int num = 1;
            for (int i = 0; i < scene.areaLightBuffer.size() && scene.materials[mtlHandle.index()].hasProperty("diffuseColor"); i++) {
                auto p1 = (rand() % 1000) / 1000.f;
                auto p2 = (rand() % 1000) / 1000.f;
                auto lightRay = Ray{ hitObject->hitPoint ,glm::normalize(scene.areaLightBuffer[i].position
                    + p1 * scene.areaLightBuffer[i].u
                    + p2 * scene.areaLightBuffer[i].v
                    - hitObject->hitPoint) };
                auto tempHitObject = closestHitObject(lightRay);
                auto [t, emitted] = closestHitLight(lightRay);
                if (tempHitObject && tempHitObject->t < t) continue;
                float A = scene.areaLightBuffer[i].A;
                auto cos_0 = abs(glm::dot(lightRay.direction, scene.areaLightBuffer[i].n));
                auto n_dot_in = abs(glm::dot(hitObject->normal, lightRay.direction));
                Lo += attenuation * n_dot_in * (emitted * cos_0 / (t * t)) / (1.f / A);
                num++;
            }
            return emitted + Lo / (1.f * num);
        }
        else if (t != FLOAT_INF) {
            return emitted;
        }
        else {
            return Vec3{0};
        }
    }
    RGB PathPhotoMappingRenderer::trace_indirect(const Ray& r, int currDepth) {
        auto hitObject = closestHitObject(r);
        auto [ t, emitted ] = closestHitLight(r);
        // hit object
        if (hitObject && hitObject->t < t) {
            auto mtlHandle = hitObject->material;
            auto scattered = shaderPrograms[mtlHandle.index()]->shade(r, hitObject->hitPoint, hitObject->normal);
            auto scatteredRay = scattered.ray;
            auto attenuation = scattered.attenuation;
            auto emitted = scattered.emitted;

            if (scene.materials[mtlHandle.index()].hasProperty("ior")) {
                RGB next = Vec3(0.f);
                if (currDepth > depth && rand() % 1000 < 200) next = (attenuation + scattered.refraction) * scene.ambient.constant;
                else {
                    auto reflex = trace_indirect(scattered.ray, currDepth + 1);
                    RGB refraction = Vec3(0.f);
                    if (scattered.refraction_ray.direction != Vec3(0.f)) refraction = trace_indirect(scattered.refraction_ray, currDepth + 1);
                    next = attenuation * reflex + refraction * scattered.refraction;
                }
                return emitted + (next / scattered.pdf)  / 0.8f;
            }

            Vec3 dir{ 0,0,0 };
            Vec3 r{ 0,0,0 };

            float max_r = 0.f;

            auto k_nearst = photon_tree.searchKNearest(photon{ hitObject->hitPoint }, 50, max_r);
            
            for (int i = 0; i < k_nearst.size(); i++) {
                auto p = photons[k_nearst[i]];
                float n_dot_in = glm::dot(-(p.in.direction), hitObject->normal);
                if (n_dot_in <= 0.f) continue;
                dir += -(p.in.direction);
                r += p.r / (float)(PI * std::pow(max_r, 2) * photonNum);
            }

            float n_dot_in = glm::dot(hitObject->normal, glm::normalize(dir));
            float pdf = scattered.pdf;

            auto Lc = (attenuation * r * n_dot_in / pdf);
            return emitted + Lc;
        }
        else if (t != FLOAT_INF) {
            return emitted;
        }
        else {
            return Vec3{0};
        }
    }

    RGB PathPhotoMappingRenderer::trace_caustics(const Ray& r, int currDepth) {
        auto hitObject = closestHitObject(r);
        auto [ t, emitted ] = closestHitLight(r);
        // hit object
        if (hitObject && hitObject->t < t) {
            auto mtlHandle = hitObject->material;
            auto scattered = shaderPrograms[mtlHandle.index()]->shade(r, hitObject->hitPoint, hitObject->normal);
            auto scatteredRay = scattered.ray;
            auto attenuation = scattered.attenuation;
            auto emitted = scattered.emitted;

            if (scene.materials[mtlHandle.index()].hasProperty("ior")) {
                RGB next = Vec3(0.f);
                if (currDepth > depth && rand() % 1000 < 200) next = (attenuation + scattered.refraction) * scene.ambient.constant;
                else {
                    auto reflex = trace_caustics(scattered.ray, currDepth + 1);
                    RGB refraction = Vec3(0.f);
                    if (scattered.refraction_ray.direction != Vec3(0.f)) refraction = trace_caustics(scattered.refraction_ray, currDepth + 1);
                    next = attenuation * reflex + refraction * scattered.refraction;
                }
                return emitted + (next / scattered.pdf)  / 0.8f;
            }

            Vec3 dir{ 0,0,0 };
            Vec3 r{ 0,0,0 };

            float max_r = 0.f;

            auto k_nearst = photon_tree.searchKNearest(photon{ hitObject->hitPoint }, 50, max_r);
            
            for (int i = 0; i < k_nearst.size(); i++) {
                auto p = photons[k_nearst[i]];
                float n_dot_in = glm::dot(-(p.in.direction), hitObject->normal);
                if (n_dot_in <= 0.f) continue;
                dir += -(p.in.direction);
                r += p.r / (float)(PI * std::pow(max_r, 2) * photonNum);
            }

            float n_dot_in = glm::dot(hitObject->normal, glm::normalize(dir));
            float pdf = scattered.pdf;

            auto Lc = (attenuation * r * n_dot_in / pdf);
            return emitted + Lc;
        }
        else if (t != FLOAT_INF) {
            return emitted;
        }
        else {
            return Vec3{0};
        }
    }
}