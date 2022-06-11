#pragma once
#ifndef __SIMPLE_PATH_TRACER_HPP__
#define __SIMPLE_PATH_TRACER_HPP__

#include "scene/Scene.hpp"
#include "Ray.hpp"
#include "Camera.hpp"
#include "intersections/HitRecord.hpp"

#include "shaders/ShaderCreator.hpp"

#include "KdTree.hpp"

#include <tuple>
#include <map>

namespace PhotoMapping
{
    using namespace NRenderer;
    using namespace std;

    class PhotoMappingRenderer
    {
    public:
    private:
        SharedScene spScene;
        Scene& scene;

        unsigned int width;
        unsigned int height;
        unsigned int depth;
        unsigned int samples;
        unsigned int photonNum;

        using SCam = PhotoMapping::Camera;
        SCam camera;

        vector<SharedShader> shaderPrograms;
    public:
        PhotoMappingRenderer(SharedScene spScene)
            : spScene               (spScene)
            , scene                 (*spScene)
            , camera                (spScene->camera)
        {
            width = scene.renderOption.width;
            height = scene.renderOption.height;
            depth = scene.renderOption.depth;
            samples = scene.renderOption.samplesPerPixel;
            photonNum = scene.renderOption.photonNum;
        }
        ~PhotoMappingRenderer() = default;

        using RenderResult = tuple<RGBA*, unsigned int, unsigned int>;
        RenderResult render();
        void release(const RenderResult& r);

    private:
        void renderTask(RGBA* pixels, int width, int height);
        void pthotoMap();
        void pthotoMap(Ray& ray, Vec3 r, int depth);
        std::vector<photon> photons;
        KdTree tree;
        KDTree<3> photon_tree;
        RGB gamma(const RGB& rgb);
        RGB trace(const Ray& ray, int currDepth, bool gather);
        HitRecord closestHitObject(const Ray& r);
        tuple<float, Vec3> closestHitLight(const Ray& r);
    };
}

#endif