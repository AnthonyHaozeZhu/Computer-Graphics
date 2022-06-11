#include "shaders/Glass.hpp"
#include "samplers/SamplerInstance.hpp"

namespace SimplePathTracer
{
    Glass::Glass(Material& material, vector<Texture>& textures)
        : Shader                (material, textures)
    {
        auto ior = material.getProperty<Property::Wrapper::FloatType>("ior");
        this->ior = (*ior).value;

        auto absorbed = material.getProperty<Property::Wrapper::RGBType>("absorbed");
        this->absorbed = (*absorbed).value;
    }
    Scattered Glass::shade(const Ray& ray, const Vec3& hitPoint, const Vec3& normal) const {
        Vec3 origin = hitPoint;
        Vec3 n = glm::normalize(normal);
        float i = ior;
        Vec3 in = glm::normalize(ray.direction);
        float in_cos = glm::dot(in, n);
        if(in_cos > 0.f) n = -n, i = 1.f / i;
        Vec3 reflex_dir = glm::normalize(in + n);

        float pdf = 1.f;

        float R0 = std::pow((i - 1.f) / (i + 1.f), 2);

        Vec3 reflex = (R0 + (1.f - R0) * (float)std::pow(1 - std::abs(in_cos), 5)) * absorbed;
        auto refraction = absorbed - reflex;

        Vec3 x = glm::normalize(reflex_dir + in);
        Vec3 y = glm::normalize(-n);

        float sin_r = std::pow(1 - std::pow(in_cos, 2), 0.5) / i;
        float cos_r = std::pow(1 - std::pow(sin_r, 2), 0.5);

        auto refraction_dir = glm::normalize(x * sin_r + y * cos_r);
        if (sin_r > 1.f) {
            reflex = absorbed; 
            refraction = Vec3(0.f);
            refraction_dir = Vec3(0.f);
        }
        return {
            Ray{origin, reflex_dir},
            reflex,
            Vec3{0},
            pdf,
            Ray{origin, refraction_dir},
            refraction
        };
    }
}