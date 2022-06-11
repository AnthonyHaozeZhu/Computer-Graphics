#include "server/Server.hpp"
#include "scene/Scene.hpp"
#include "component/RenderComponent.hpp"
#include "Camera.hpp"

#include "PhotonMapping.hpp"

using namespace std;
using namespace NRenderer;

namespace PhotoMapping
{
    class Adapter : public RenderComponent
    {
        void render(SharedScene spScene) {
            PhotoMappingRenderer renderer{spScene};
            auto renderResult = renderer.render();
            auto [ pixels, width, height ]  = renderResult;
            getServer().screen.set(pixels, width, height);
            renderer.release(renderResult);
        }
    };
}

const static string description = 
    "photo mapping\n";

REGISTER_RENDERER(PhotoMapping, description, PhotoMapping::Adapter);