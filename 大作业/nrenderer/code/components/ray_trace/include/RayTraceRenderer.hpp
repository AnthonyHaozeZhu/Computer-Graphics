#pragma once
#ifndef __RAY_CAST_HPP__
#define __RAY_CAST_HPP__

#include "scene/Scene.hpp"

#include "Camera.hpp"
#include "intersections/intersections.hpp"

#include "shaders/ShaderCreator.hpp"

#include "KdTree.hpp"

namespace RayTrace
{
    using namespace NRenderer;
    class RayTraceRenderer
    {
    private:
        SharedScene spScene;
        Scene& scene;
        RayTrace::Camera camera;

        vector<SharedShader> shaderPrograms;
    public:
        RayTraceRenderer(SharedScene spScene)
            : spScene               (spScene)
            , scene                 (*spScene)
            , camera                (spScene->camera)
        {}
        ~RayTraceRenderer() = default;

        using RenderResult = tuple<RGBA*, unsigned int, unsigned int>;
        RenderResult render();
        void release(const RenderResult& r);
        KDTree tree;
    private:
        RGB gamma(const RGB& rgb);
        RGB trace(const Ray& r, int depth);
        HitRecord closestHit(const Ray& r);
    };
}

#endif