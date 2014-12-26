#include <QColor>
#include <QApplication>

#include <cg/visualization/viewer_adapter.h>
#include <cg/visualization/draw_util.h>

#include <cg/io/point.h>

#include <cg/dcel/kirkpatrick.h>

#include <queue>
#include <iostream>
#include <cmath>

using cg::DCEL;
using cg::kirkpatrick_localization;
using cg::line;
using cg::vertex;
using cg::edge;
using cg::triangle_k;
using cg::node_triangle;

using cg::point_2f;
using cg::point_2;

struct dcel_viewer : cg::visualization::viewer_adapter
{
    dcel_viewer()
        : kl(line{0, -1, 0}, line{1, 0, 0})
        , points_seted(0)
        , mode(line_mode)
        , in_color(false), show_triangulation(false)
        , lvl(0)
    {
        kl.build_triangulation(deleted_vertices);
    }

    void draw_shift_edge(const point_2 &  v, const point_2 & u, QColor edge_color,
                         cg::visualization::drawer_type & drawer) const
    {
        double a = u.y - v.y, b = v.x - u.x,
               nx = a / std::sqrt(a * a + b * b),
               ny = b / std::sqrt(a * a + b * b);

        double d = std::sqrt((v.x - u.x) * (v.x - u.x) + (v.y - u.y) * (v.y - u.y));
        double cosf = (u.x - v.x) / d, sinf = (u.y - v.y) / d;

        double r = 10, nr = 2;
        auto p1 = point_2{v.x + cosf * r - nx * nr, v.y + sinf * r - ny * nr};
        auto p2 = point_2{u.x - cosf * r - nx * nr, u.y - sinf * r - ny * nr};

        drawer.set_color(edge_color);
        drawer.draw_line(p1, p2);
        drawer.set_color(Qt::white);
        drawer.draw_point(p2, 4);
    }

    void draw_ray(const std::shared_ptr<edge> & e, QColor ray_color,
                  cg::visualization::drawer_type & drawer) const
    {
        point_2 v;
        if (e->origin->line1) {
            v = cg::intersection_point<double>(*e->origin->line1, *e->origin->line2);
        } else {
            v = cg::intersection_point<double>(*e->next->origin->line1, *e->next->origin->line2);
        }
        double x = -e->line_link->b, y = e->line_link->a;
        double px = x / std::sqrt(x * x + y * y);
        double py = y / std::sqrt(x * x + y * y);

        double r = 500;
        point_2 u = point_2{v.x + r * px, v.y + r * py};

        if (e->origin->line1) {
            draw_shift_edge(v, u, ray_color, drawer);
        } else {
            draw_shift_edge(u, v, ray_color, drawer);
        }
    }

    void draw_red_ray(const point_2 & p, const point_2 & s, float width,
                      cg::visualization::drawer_type & drawer) const
    {
        double sx = s.x / std::sqrt(s.x * s.x + s.y * s.y);
        double sy = s.y / std::sqrt(s.x * s.x + s.y * s.y);

        double r = 500;
        point_2 q = point_2{p.x + r * sx, p.y + r * sy};

        drawer.draw_line(p, q, width);
    }

    inline line line_by_points(const point_2 & v, const point_2 & u) const
    {
        double a = u.y - v.y, b = v.x - u.x, c = -(a * v.x + b * v.y);
        line l(a, b, c);
        if (!cg::is_direct_vector_right(l)) {
            l.inverse_vector();
        }
        return l;
    }

    void bfs_draw(cg::visualization::drawer_type & drawer, const DCEL & dcel) const
    {
        auto visited = std::vector<int>(dcel.max_vertex);
        std::queue<std::shared_ptr<vertex>> q;
        auto s = dcel.inf_node;
        visited[s->id] = 1;
        q.push(s);

        while (!q.empty()) {
            auto v = q.front();
            q.pop();

            if (v->line1) { // not inf_node
                drawer.set_color(Qt::red);
                drawer.draw_point(cg::intersection_point<double>(*v->line1, *v->line2), 5);
            }

            auto out_edge = v->e;
            do {
                //std::cout << out_edge->id << std::endl;
                auto u = out_edge->next->origin;
                if (!visited[u->id]) {
                    q.push(u);
                    visited[u->id] = 1;
                }
                if (!v->line1 || !u->line1) {
                    draw_ray(out_edge, Qt::yellow, drawer);
                } else {
                    auto p1 = cg::intersection_point<double>(*v->line1, *v->line2);
                    auto p2 = cg::intersection_point<double>(*u->line1, *u->line2);
                    if (out_edge->hull_edge) {
                        draw_shift_edge(p1, p2, Qt::darkGreen, drawer);
                    } else if (out_edge->triangle_edge) {
                        drawer.set_color(QColor(255, 69, 0));
                        drawer.draw_line(p1, p2, 1);
                    } else {
                        draw_shift_edge(p1, p2, Qt::yellow, drawer);
                    }
                }
                out_edge = out_edge->twin->next;
            } while (out_edge->id != v->e->id);
        }
    }

    void draw_blue_line(cg::visualization::drawer_type & drawer) const
    {
        double r = 1000;
        double d = std::sqrt((l1.x - l2.x) * (l1.x - l2.x) + (l1.y - l2.y) * (l1.y - l2.y));
        double cosf = (l2.x - l1.x) / d;
        double sinf = (l2.y - l1.y) / d;

        auto p1 = point_2{l1.x - cosf * r, l1.y - sinf * r};
        auto p2 = point_2{l2.x + cosf * r, l2.y + sinf * r};

        drawer.set_color(QColor(135, 206, 250)); // light sky blue color
        drawer.draw_line(l1, l2);
        drawer.draw_line(p1, l1);
        drawer.draw_line(l2, p2);
    }

    void draw_triangle(const triangle_k & t, cg::visualization::drawer_type & drawer) const
    {
        drawer.set_color(QColor(175, 238, 238)); // paleturquoise
        for (int i = 0; i < 3; i++) {
            point_2 p1 = cg::intersection_point<double>(*t[i].l1, *t[i].l2);
            point_2 p2 = cg::intersection_point<double>(*t[(i + 1) % 3].l1, *t[(i + 1) % 3].l2);
            drawer.draw_line(p1, p2, 3.3);
        }
    }

    void draw_triangle_level(cg::visualization::drawer_type & drawer) const
    {
        int rem = kl.levels.size() - lvl - 1;
        line l1(1, 0, -face_point->x), l2(0, 1, -face_point->y);
        std::shared_ptr<node_triangle> node = kl.root;

        while (!node->is_leaf && node->depth > kl.root->depth - rem) {
            for (auto child : node->children) {
                if (cg::triangle_contains_convex_point(*child->t, l1, l2)) {
                    node = child;
                    break;
                }
            }
        }
        draw_triangle(*node->t, drawer);
    }

    void draw(cg::visualization::drawer_type & drawer) const
    {
        if (show_triangulation) {
            bfs_draw(drawer, kl.levels[lvl].dcel);
        } else {
            bfs_draw(drawer, kl.dcel);
        }

        if (face_point) {
            drawer.set_color(QColor(173, 255, 47)); // greenyellow
            drawer.draw_point(*face_point, 3);
        }

        if (points_seted == 2) {
            draw_blue_line(drawer);
            if (in_color) {
                std::vector<std::shared_ptr<edge>> intersected_edges;
                kl.dcel.get_intersected_edges(line_by_points(l1, l2), intersected_edges);

                drawer.set_color(Qt::blue);
                for (auto e : intersected_edges) {
                    auto vu = e;
                    do {
                        auto v = vu->origin;
                        auto u = vu->next->origin;
                        if (!v->line1 || !u->line1) {
                            draw_ray(vu, Qt::blue, drawer);
                        } else {
                            auto p1 = cg::intersection_point<double>(*v->line1, *v->line2);
                            auto p2 = cg::intersection_point<double>(*u->line1, *u->line2);
                            draw_shift_edge(p1, p2, Qt::blue, drawer);
                        }
                        vu = vu->next;
                    } while (vu->id != e->id);
                }
            }
        }

        if (face_edge) {
            if (!show_triangulation) {
                auto e = face_edge;
                do {
                    if (!e->origin->line1 || !e->next->origin->line1) {
                        draw_ray(e, Qt::blue, drawer);
                    } else {
                        auto p1 = cg::intersection_point<double>(*e->origin->line1, *e->origin->line2);
                        auto p2 = cg::intersection_point<double>(*e->next->origin->line1, *e->next->origin->line2);
                        draw_shift_edge(p1, p2, Qt::blue, drawer);
                    }
                    e = e->next;
                } while (e->id != face_edge->id);
            } else {
                if (face_point) {
                    draw_triangle_level(drawer);
                }
            }
        }

        if (show_triangulation) {
            for (auto v : deleted_vertices[lvl]) {
                auto p = cg::intersection_point<double>(*v->line1, *v->line2);
                drawer.set_color(Qt::darkMagenta);
                drawer.draw_point(p, 10);
            }
        }
    }

    void print(cg::visualization::printer_type & p) const
    {
        p.corner_stream() << "press 'l' for 'Line' mode, 'f' for 'Face', 'c' or 'h' to clear, 'b' for blue color, 't' for triangles" << cg::visualization::endl
                          << "in Line mode: press mouse rbutton to add point then press one more time to draw line" << cg::visualization::endl
                          << "current mode: ";
        switch (mode)
        {
            case line_mode : p.corner_stream() << "Line" << cg::visualization::endl; break;
            case face_mode : p.corner_stream() << "Face" << cg::visualization::endl; break;
        }

        if (in_color) {
            p.corner_stream() << "Color: ON" << cg::visualization::endl;
        } else {
            p.corner_stream() << "Color: OFF" << cg::visualization::endl;
        }

        if (show_triangulation) {
            p.corner_stream() << "Triangles: ON" << cg::visualization::endl;
        } else {
            p.corner_stream() << "Triangles: OFF" << cg::visualization::endl;
        }

        p.corner_stream() << "Level: " << (lvl + 1) << " (" << kl.levels.size() << ")" << cg::visualization::endl;

        if (face_point) {
            p.corner_stream() << "Point: " << (*face_point) << cg::visualization::endl;
        }
    }

    bool on_release(const point_2f & p)
    {
        if (mode == line_mode) {
            switch (points_seted) {
            case 0:
                points_seted = 1;
                l1 = l2 = p;
                break;
            case 1:
                points_seted = 0;
                break;
            case 2:
                l2 = p;
                kl.add_line(line_by_points(l1, l2));
                kl.build_triangulation(deleted_vertices);
                points_seted = 0;
                break;
            }
        } else { // mode == face_mode
            face_point = std::make_shared<point_2>(p);
            if (!show_triangulation) {
                face_edge = kl.naive_localization(p);
            } else {
                face_edge = kl.fast_localization(p);
            }
        }

        return true;
    }

    bool on_move(const point_2f & p)
    {
        if (points_seted > 0) {
            l2 = p;
            double current_r = (l1.x - l2.x) * (l1.x - l2.x) + (l1.y - l2.y) * (l1.y - l2.y);

            if (current_r > 25) {
                points_seted = 2;
            } else {
                points_seted = 1;
            }
        }

        return points_seted > 0;
    }


    bool on_key(int key)
    {
        switch (key) {
        case Qt::Key_C:
            kl = kirkpatrick_localization(line{1, -1, 0}, line{-1, -1, 0});
            kl.build_triangulation(deleted_vertices);
            points_seted = 0;
            face_edge = std::shared_ptr<edge>();
            face_point = nullptr;
            break;
        case Qt::Key_H:
            kl = kirkpatrick_localization(line{0, -1, 0}, line{1, 0, 0});
            kl.build_triangulation(deleted_vertices);
            points_seted = 0;
            face_edge = std::shared_ptr<edge>();
            face_point = nullptr;
            break;
        case Qt::Key_T:
            show_triangulation = !show_triangulation;
            face_point = nullptr;
            break;
        case Qt::Key_L:
            mode = line_mode;
            face_edge = std::shared_ptr<edge>();
            face_point = nullptr;
            break;
        case Qt::Key_F:
            mode = face_mode;
            points_seted = 0;
            face_point = nullptr;
            break;
        case Qt::Key_B:
            in_color = !in_color;
            break;
        case Qt::Key_M:
            kl = kirkpatrick_localization(kl);
            kl.build_triangulation(deleted_vertices);
            break;
        case Qt::Key_Right:
            lvl = std::min(lvl + 1, int(kl.levels.size()) - 1);
            break;
        case Qt::Key_Left:
            lvl = std::max(0, lvl - 1);
            break;
        default:
            return false;
        }

        return true;
    }

private:
    kirkpatrick_localization kl;

    int points_seted;
    point_2 l1, l2; // line points
    std::shared_ptr<point_2> face_point;

    enum point_mode
    {
        line_mode = 0,
        face_mode = 1
    };

    point_mode mode;
    std::shared_ptr<edge> face_edge;

    bool in_color;
    bool show_triangulation;

    int lvl;
    std::vector< std::vector< std::shared_ptr<vertex> > > deleted_vertices;
};

int main(int argc, char ** argv)
{
    QApplication app(argc, argv);
    dcel_viewer viewer;
    cg::visualization::run_viewer(&viewer, "Doubly-Connected Edge List");
}
