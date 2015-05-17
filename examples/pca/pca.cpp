#include <QColor>
#include <QApplication>
#include <QTextStream>

#include <cg/visualization/viewer_adapter.h>
#include <cg/visualization/draw_util.h>

#include <cg/pca/pca.h>
#include <cg/pca/weights.h>

#include <cg/io/point.h>

#include <boost/optional.hpp>

#include <iostream>
#include <vector>
#include <set>

using cg::point_2f;
using cg::point_2;
using cg::weighted_pca;
using cg::weight_type;

struct pca_viewer : cg::visualization::viewer_adapter
{
    pca_viewer()
        : pca_type(cg::standard_pca),
          eig_ret(boost::none),
          drag_cursor(boost::none)
    {
        wpca.set_weight(cg::unit_weight());
    }

    void draw(cg::visualization::drawer_type & drawer) const
    {
        drawer.set_color(Qt::red);
        drawer.draw_line(point_2{0, -300}, point_2{0, 300});
        drawer.draw_line(point_2{-300, 0}, point_2{300, 0});

        drawer.set_color(Qt::green);
        for (point_2 p : points) {
            drawer.draw_point(p, 2);
        }

        // draw point for remove
        if (cur_cursor != boost::none) {
            drawer.set_color(Qt::yellow);
            drawer.draw_point(cur_cursor.get(), 5);
        }

        if (eig_ret != boost::none) {
            double r = 200;
            std::pair<point_2, point_2> eig = eig_ret.get();
            point_2 m = eig.second, d = eig.first;
            point_2 from = point_2{m.x - r * d.x, m.y - r * d.y};
            point_2 to   = point_2{m.x + r * d.x, m.y + r * d.y};

            drawer.set_color(Qt::blue);
            drawer.draw_line(from, to, 2);
        }
    }

    void print(cg::visualization::printer_type & p) const
    {
        p.corner_stream() << "pca visualization" << cg::visualization::endl
                          << "press 's' for PCA (default)" << cg::visualization::endl
                          << "press 'h' for manhattan WPCA" << cg::visualization::endl
                          << "press 'e' for euclidian WPCA" << cg::visualization::endl
                          << "press 'm' for mahalanobis WPCA" << cg::visualization::endl;

        p.corner_stream() << "current mode: ";
        switch (pca_type) {
            case cg::standard_pca: p.corner_stream() << "standard"; break;
            case cg::manhattan_pca: p.corner_stream() << "manhattan"; break;
            case cg::euclidian_pca: p.corner_stream() << "euclidian"; break;
            case cg::mahalanobis_pca: p.corner_stream() << "mahalanobis"; break;
        }
        p.corner_stream() << cg::visualization::endl;
    }

    bool on_press(const point_2f & p)
    {
        if (cur_cursor != boost::none) {
            drag_cursor = cur_cursor.get();
            beg_cursor = cur_cursor.get();
            return true;
        }

        return false;
    }

    bool on_release(const point_2f & p)
    {
       if (cur_cursor == boost::none) {
           points.insert(p);
           cur_cursor = p;
       } else if (beg_cursor == cur_cursor) {
           points.erase(cur_cursor.get());
           cur_cursor = boost::none;
       }
       drag_cursor = boost::none;

       if (points.size() > 1) {
           eig_ret = wpca.find_principal_component(points);
       }

       return true;
    }

    bool on_double_click(const point_2f & p)
    {
       points.clear();
       cur_cursor = boost::none;
       eig_ret = boost::none;

       return true;
    }

    bool on_key(int key)
    {
        switch (key) {
            case Qt::Key_S:
                pca_type = cg::standard_pca;
                wpca.set_weight(cg::unit_weight());
                break;
            case Qt::Key_H:
                pca_type = cg::manhattan_pca;
                wpca.set_weight(cg::manhattan_weight());
                break;
            case Qt::Key_E:
                pca_type = cg::euclidian_pca;
                wpca.set_weight(cg::euclidian_weight());
                break;
            case Qt::Key_M:
                pca_type = cg::mahalanobis_pca;
                wpca.set_weight(cg::mahalanobis_weight(points));
                break;
            default: return false;
        }

        if (points.size() > 1) {
            eig_ret = wpca.find_principal_component(points);
        }
        return true;
    }

    void set_cur_cursor(const point_2 & pnt)
    {
        cur_cursor = boost::none;
        double max_r = 0;

        for (point_2 p : points) {
            double current_r = (p.x - pnt.x) * (p.x - pnt.x) + (p.y - pnt.y) * (p.y - pnt.y);

            if ((cur_cursor == boost::none && current_r < 100) || (cur_cursor != boost::none && current_r < max_r)) {
                cur_cursor = p;
                max_r = current_r;
            }
        }

//        return cur_cursor != boost::none;
    }

    bool on_move(const point_2f & p)
    {
        if (drag_cursor != boost::none) {
            points.erase(drag_cursor.get());
            drag_cursor = p;
            points.insert(drag_cursor.get());
            if (points.size() > 1) {
                eig_ret = wpca.find_principal_component(points);
            }
        }
        set_cur_cursor(p);
        return true;
    }

private:
    std::set<point_2> points; // all points on the plane
    boost::optional<point_2> cur_cursor; // current cursor (for remove)
    boost::optional<point_2> drag_cursor;
    boost::optional<point_2> beg_cursor;

    weight_type pca_type;
    weighted_pca wpca;
    boost::optional<std::pair<point_2, point_2>> eig_ret;
};

int main(int argc, char ** argv)
{
   QApplication app(argc, argv);
   pca_viewer viewer;
   cg::visualization::run_viewer(&viewer, "PCA");
}
