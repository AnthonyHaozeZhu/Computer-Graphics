#include "server/Server.hpp"

#include "PhotonMapping.hpp"

#include "VertexTransformer.hpp"
#include "intersections/intersections.hpp"

#include "glm/gtc/matrix_transform.hpp"
#include "Onb.hpp"
#include <omp.h>

#define KDTREE

namespace PhotoMapping
{
    RGB PhotoMappingRenderer::gamma(const RGB& rgb) {
        return glm::sqrt(rgb);
    }
    int index = 0;

    Vec3 randomDir(const Vec3& normal) {
        Vec3 random = defaultSamplerInstance<HemiSphere>().sample3d();
        Onb onb{ normal };
        Vec3 direction = glm::normalize(onb.local(random));
        return direction;
    }

    void PhotoMappingRenderer::pthotoMap() {
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
        std::vector<VecXf<3>> points(photons.size());
        for (int i = 0; i < photons.size(); i++) {
            for (int j = 0; j < 3; j++) {
                points[i][j] = photons[i].position[j];
            }
        }
        photon_tree.buildTree(std::move(points));
        sprintf(log, "photon num: %u", photons.size());
        getServer().logger.log(log);
    }

    void PhotoMappingRenderer::pthotoMap(Ray& ray, Vec3 r, int currDepth) {
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

    void PhotoMappingRenderer::renderTask(RGBA* pixels, int width, int height) {
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
                    color += trace(ray, 0, 1);
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

    auto PhotoMappingRenderer::render() -> RenderResult {
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
        PhotoMappingRenderer::renderTask(pixels, width, height);
        getServer().logger.log("Done...");
        return {pixels, width, height};
    }

    void PhotoMappingRenderer::release(const RenderResult& r) {
        auto [p, w, h] = r;
        delete[] p;
    }

    HitRecord PhotoMappingRenderer::closestHitObject(const Ray& r) {
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
    
    tuple<float, Vec3> PhotoMappingRenderer::closestHitLight(const Ray& r) {
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

    RGB PhotoMappingRenderer::trace(const Ray& r, int currDepth, bool gather) {
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
                    auto reflex = trace(scattered.ray, currDepth + 1, gather);
                    RGB refraction = Vec3(0.f);
                    if (scattered.refraction_ray.direction != Vec3(0.f)) refraction = trace(scattered.refraction_ray, currDepth + 1, gather);
                    next = attenuation * reflex + refraction * scattered.refraction;
                }
                return emitted + (next / scattered.pdf)  / 0.8f;
            }

            Vec3 dir{ 0,0,0 };
            Vec3 r{ 0,0,0 };

            float max_r = 0.f;
            auto k_nearst = photon_tree.nearstSearch(hitObject->hitPoint, max_r, 30);

            
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