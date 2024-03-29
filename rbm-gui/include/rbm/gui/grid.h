#include "rbm/gui/entity.h"
#include "rbm/gui/grid_renderer.h"

namespace rbm {

class GridEntity: public Entity {
public:
    GridEntity():
        grid(-1)
    {}
    bool init(Renderers& renderers) override
    {
        auto renderer = renderers.get<GridRenderer>();
        grid = renderer->create_grid();
        renderer->configure_grid(grid, 10, 0.03, 10, 10);
        return true;
    }
    void render(Renderers& renderers) override
    {
        auto renderer = renderers.get<GridRenderer>();
        renderer->queue_grid(grid);
    }

private:
    int grid;
};

} // namespace rbm
