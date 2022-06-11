#include "server/Server.hpp"

#include "SimplePathTracer.hpp"

#include "VertexTransformer.hpp"
#include "intersections/intersections.hpp"

#include "glm/gtc/matrix_transform.hpp"

#include <omp.h>

#define KDTREE

namespace SimplePathTracer
{
    RGB SimplePathTracerRenderer::gamma(const RGB& rgb) {
        return glm::sqrt(rgb);
    }

    constexpr auto E = 2.718281828459f;

    static inline Vec3 tinyChange(Vec3& dir) {
        Vec3 d = dir;
        float s1 = 1.f / 1024.f;
        float s2 = 1.f / 16.f;
        for (int i = 0; i < 3; i++) {
            float dx = s2 * std::pow(E, -std::log(s2 / s1)) * ((rand() % 10000 - 5000) / 10000.f);
            d[i] += dx;
        }
        dir = glm::normalize(d);
        return dir;
    }

    static inline Ray randomRay(int height, int width, int i, int j, auto& camera) {
        auto r = defaultSamplerInstance<UniformInSquare>().sample2d();
        float rx = r.x;
        float ry = r.y;
        float x = (float(j) + rx) / float(width);
        float y = (float(i) + ry) / float(height);
        auto ray = camera.shoot(x, y);
        return ray;
    }

    static inline float calIuminance(const RGB& color) {
        return 0.299f * color[0] + 0.587f * color[1] * 0.114f * color[2];
    }

    void SimplePathTracerRenderer::renderTask(RGBA* pixels, int width, int height) {
        #pragma omp parallel for
        for(int i=0; i<height; i++) {
            for (int j=0; j<width; j++) {
                Vec3 color{0, 0, 0};
                for (int k = 0; k < samples; k++) {
                    auto ray = randomRay(height, width, i, j, camera);
                    color += trace(ray, 0);
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

    auto SimplePathTracerRenderer::render() -> RenderResult {
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

        SimplePathTracerRenderer::renderTask(pixels, width, height);
        getServer().logger.log("Done...");
        return {pixels, width, height};
    }

    void SimplePathTracerRenderer::release(const RenderResult& r) {
        auto [p, w, h] = r;
        delete[] p;
    }

    HitRecord SimplePathTracerRenderer::closestHitObject(const Ray& r) {
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
    
    tuple<float, Vec3> SimplePathTracerRenderer::closestHitLight(const Ray& r) {
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

    RGB SimplePathTracerRenderer::trace(const Ray& r, int currDepth) {
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
            if(scene.materials[mtlHandle.index()].hasProperty("diffuseColor")){
                float n_dot_in = glm::dot(hitObject->normal, scatteredRay.direction);
                if (currDepth > depth && rand() % 1000 < 200) next = attenuation * n_dot_in * scene.ambient.constant;
                else next = attenuation * n_dot_in * trace(scatteredRay, currDepth + 1);
                
            }
            else if (scene.materials[mtlHandle.index()].hasProperty("ior")) {
                if (currDepth > depth && rand() % 1000 < 200) next = (attenuation + scattered.refraction) * scene.ambient.constant;
                else {
                    auto reflex = trace(scatteredRay, currDepth + 1);
                    RGB refraction = Vec3(0.f);
                    if (scattered.refraction_ray.direction != Vec3(0.f))  refraction = trace(scattered.refraction_ray, currDepth + 1);
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
}