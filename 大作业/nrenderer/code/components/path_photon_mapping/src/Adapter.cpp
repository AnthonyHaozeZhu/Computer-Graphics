#include "server/Server.hpp"
#include "scene/Scene.hpp"
#include "component/RenderComponent.hpp"
#include "Camera.hpp"

#include "PathPhotonMapping.hpp"

using namespace std;
using namespace NRenderer;

namespace PathPhotonMapping
{
    class Adapter : public RenderComponent
    {
        void render(SharedScene spScene) {
            PathPhotoMappingRenderer renderer{spScene};
            auto renderResult = renderer.render();
            auto [ pixels, width, height ]  = renderResult;
            getServer().screen.set(pixels, width, height);
            renderer.release(renderResult);
        }
    };
}

const static string description = 
    "photo mapping\n";

REGISTER_RENDERER(PathPhotonMapping, description, PathPhotonMapping::Adapter);