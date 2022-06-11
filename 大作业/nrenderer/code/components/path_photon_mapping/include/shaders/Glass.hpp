#pragma once
#ifndef __GLASS_HPP__
#define __GLASS_HPP__

#include "Shader.hpp"

namespace PathPhotonMapping
{
    class Glass : public Shader
    {
    private:
        Vec3 absorbed;
        float ior;
    public:
        Glass(Material& material, vector<Texture>& textures);
        Scattered shade(const Ray& ray, const Vec3& hitPoint, const Vec3& normal) const;
    };
}

#endif