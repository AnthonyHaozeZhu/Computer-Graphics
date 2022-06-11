#include "RayTraceRenderer.hpp"

#include "VertexTransformer.hpp"
#include "intersections/intersections.hpp"

#include "omp.h"

#define KDTREE

namespace RayTrace
{
    void RayTraceRenderer::release(const RenderResult& r) {
        auto [p, w, h] = r;
        delete[] p;
    }
    RGB RayTraceRenderer::gamma(const RGB& rgb) {
        return glm::sqrt(rgb);
    }
    auto RayTraceRenderer::render() -> RenderResult {
        auto width = scene.renderOption.width;
        auto height = scene.renderOption.height;
        auto pixels = new RGBA[width*height];

        VertexTransformer vertexTransformer{};
        vertexTransformer.exec(spScene);
        
        ShaderCreator shaderCreator{};
        for (auto& mtl : scene.materials) {
            shaderPrograms.push_back(shaderCreator.create(mtl, scene.textures));
        }
        tree.Insert(scene.sphereBuffer, scene.triangleBuffer);
        #pragma omp parallel for
        for (int i=0; i<height; i++) {
            for (int j=0; j < width; j++) {
                auto ray = camera.shoot(float(j) / float(width), float(i) / float(height));
                auto color = trace(ray, 0);
                color = clamp(color);
                color = gamma(color);
                pixels[(height - i - 1) * width + j] = { color, 1 };
            }
        }

        return {pixels, width, height};
    }
    
    RGB RayTraceRenderer::trace(const Ray& r, int depth) {
        if ((depth >= 20 && rand() % 100 < 20) || scene.pointLightBuffer.size() < 1) return {0, 0, 0};
        auto& l = scene.pointLightBuffer[0];
        auto closestHitObj = closestHit(r);
        if (closestHitObj) {
            auto& hitRec = *closestHitObj;
            auto out = glm::normalize(l.position - hitRec.hitPoint);
            
            auto distance = glm::length(l.position - hitRec.hitPoint);
            auto shadowRay = Ray{hitRec.hitPoint, out};
            auto shadowHit = closestHit(shadowRay);
            auto c = shaderPrograms[hitRec.material.index()]->shade(-r.direction, out, hitRec.normal);
            auto normal = hitRec.normal;
            auto in = r.direction;
            float angle = abs(glm::dot(normal, in));
            normal = 2 * angle * normal;
            auto newRay = Ray{hitRec.hitPoint, glm::normalize(in + normal)};
            bool reflectivity = scene.materials[hitRec.material.index()].hasProperty("specularEx");
            auto c1 = shaderPrograms[hitRec.material.index()]->shade(-r.direction, glm::normalize(in + normal), hitRec.normal);
            if (glm::dot(out, hitRec.normal) < 0) {
                if (!reflectivity) return Vec3{ 0 };
                return Vec3{ 0 } + trace(newRay, depth + 1) * c1;
            }
            else if ((!shadowHit) || (shadowHit && shadowHit->t > distance)) {
                if(!reflectivity) return c * l.intensity;
                return c * l.intensity + trace(newRay, depth + 1) * c1;
            }
            else {
                if (!reflectivity) return Vec3{ 0 };
                return Vec3{ 0 } + trace(newRay, depth + 1) * c1;
            }
        }
        else {
            return Vec3{ 0 };
        }
    }

    HitRecord RayTraceRenderer::closestHit(const Ray& r) {
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
            auto hitRecord = Intersection::xPlane(r, p, 0.01, closest);
            if (hitRecord && hitRecord->t < closest) {
                closest = hitRecord->t;
                closestHit = hitRecord;
            }
        }
        return closestHit; 
    }
}