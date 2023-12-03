#include "sviz/render/entity.h"
#include "sviz/render/grid_renderer.h"

namespace sviz {

class Grid: public Entity {
public:
    Grid():
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

} // namespace sviz
