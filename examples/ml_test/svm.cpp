#include <QColor>
#include <QApplication>
#include <QTextStream>

#include <cg/visualization/viewer_adapter.h>
#include <cg/visualization/draw_util.h>

#include <cg/io/point.h>

#include <iostream>
#include <fstream>
#include <vector>

using cg::point_2f;
using cg::point_2;
using cg::segment_2;

struct quad_tree_viewer : cg::visualization::viewer_adapter
{
    quad_tree_viewer()
        : k(2.050357678467475), b(-2.5217310804908286)
    {
        std::ifstream in("LinearDataset.txt");
        if (in) {
            for (int i = 0; i < 200; i++) {

                double x, y;
                int cl;
                in >> x >> y >> cl;
                if (cl == 0) cl = -1;
                pts.push_back(point_2{x, y});
                classes.push_back(cl);
            }
        } else {
            std::cout << "oooops :(" << std::endl;
        }
    }


    void draw(cg::visualization::drawer_type & drawer) const
    {
        for (int i = 0; i < pts.size(); i++) {
            if (classes[i] == 1) {
                drawer.set_color(Qt::darkBlue);
                drawer.draw_point(pts[i], 3);
            } else {
                drawer.set_color(Qt::yellow);
                drawer.draw_point(pts[i], 3);
            }
        }

        double x1 = -20, x2 = 20;
        double y1 = k * x1 + b, y2 = k * x2 + b;

        drawer.set_color(Qt::white);
        drawer.draw_line(point_2{x1, y1}, point_2{x2, y2});

    }

    void print(cg::visualization::printer_type & p) const
    {
        p.corner_stream() << "y = " << (float) k << " * x + " << (float) b << cg::visualization::endl;
    }


private:
    std::vector<point_2> pts;
    std::vector<int> classes;

    double k, b; // k(2.1704625637247905), b(-2.3293945854743914)
};

int main(int argc, char ** argv)
{
   QApplication app(argc, argv);
   quad_tree_viewer viewer;
   cg::visualization::run_viewer(&viewer, "Points");
}
