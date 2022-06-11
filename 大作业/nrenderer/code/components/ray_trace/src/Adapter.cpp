#include "server/Server.hpp"
#include "component/RenderComponent.hpp"

#include "RayTraceRenderer.hpp"

using namespace std;
using namespace NRenderer;

namespace RayTrace
{
    class Adapter : public RenderComponent
    {
    public:
        void render(SharedScene spScene) {
            RayTraceRenderer rayTrace{spScene};
            auto result = rayTrace.render();
            auto [ pixels, width, height ] = result;
            getServer().screen.set(pixels, width, height);
            rayTrace.release(result);
        }
    };
}

const static string description = 
    "Ray Trace Renderer.\n"
    "Supported:\n"
    " - Lambertian and Phong\n"
    " - One Point Light\n"
    " - Triangle, Sphere, Plane\n"
    " - Simple Pinhole Camera\n\n"
    "Please use ray_cast.scn"
    ;

REGISTER_RENDERER(RayTrace, description, RayTrace::Adapter);
