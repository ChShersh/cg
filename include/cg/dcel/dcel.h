#pragma once

#include <cg/primitives/line_triangle.h>
#include <cg/primitives/point3d.h>
#include <cg/io/point.h>

#include <iostream>
#include <queue>

namespace cg {

    struct vertex;
    struct node_triangle;

    struct edge
    {
        int id;
        std::shared_ptr<vertex> origin;
        std::shared_ptr<edge> twin, prev, next;
        std::shared_ptr<line> line_link;
        std::shared_ptr<node_triangle> triangle_link;
        bool hull_edge;
        bool triangle_edge;

        edge()
            : hull_edge(false), triangle_edge(false)
        {}

        edge(int & num)
            : id(num), hull_edge(false), triangle_edge(false)
        {
            num++;
        }
    };

    struct vertex
    {
        int id;
        std::shared_ptr<line> line1, line2;
        std::shared_ptr<edge> e;

        vertex() {}

        vertex(int & num)
            : id(num)
        {
            num++;
        }

        vertex(const std::shared_ptr<line> & line1, const std::shared_ptr<line> & line2, int & num)
            : id(num), line1(line1), line2(line2)
        {
            num++;
        }
    };

    struct node_triangle
    {
        std::shared_ptr<triangle_k> t;
        std::shared_ptr<edge> node_edge;
        std::vector<std::shared_ptr<node_triangle>> children;
        bool is_leaf;
        int depth;

        node_triangle() {}

        node_triangle(const triangle_k & tri)
            : is_leaf(false)
        {
            t = std::make_shared<triangle_k>(tri);
        }

        node_triangle(const std::shared_ptr<triangle_k> & t)
            : t(t), is_leaf(false)
        {}
    };

    inline bool is_ray(const std::shared_ptr<edge> & e)
    {
        return !e->origin->line1 || !e->next->origin->line2;
    }

    inline std::shared_ptr<vertex> not_inf_vertex(const std::shared_ptr<edge> & e)
    {
        return (e->origin->line1 ? e->origin : e->next->origin);
    }

    struct DCEL
    {
        std::shared_ptr<vertex> inf_node;

        int max_vertex, max_edges;
        std::vector<std::shared_ptr<line>> all_lines;

        DCEL() {}

        DCEL(const line & line1, const line & line2)
            : max_vertex(0), max_edges(0)
        {
            auto l1 = std::make_shared<line>(line1);
            auto l2 = std::make_shared<line>(line2);
            if (orientation(point_2{0, 0}, point_2{-l1->b, l1->a}, point_2{-l2->b, l2->a}) == CG_RIGHT) {
                std::swap(l1, l2);
            }

            inf_node = std::make_shared<vertex>(
                std::shared_ptr<line>(), std::shared_ptr<line>(), max_vertex
            );

            auto inner_vertex = std::make_shared<vertex>(l1, l2, max_vertex);

            auto edges = std::vector<std::shared_ptr<edge>>();
            for (int i = 0; i < 8; i++) {
                edges.push_back(std::make_shared<edge>(max_edges));
            }

            for (int i = 0; i < 8; i += 2) {
                edges[i]->origin = inf_node;
                edges[i + 1]->origin = inner_vertex;
                edges[i]->twin = edges[i + 1];
                edges[i + 1]->twin = edges[i];
                edges[i]->next = edges[(i + 7) % 8];
                edges[i]->prev = edges[(i + 7) % 8];
                edges[i + 1]->next = edges[(i + 2) % 8];
                edges[i + 1]->prev = edges[(i + 2) % 8];

                if (i % 4 == 0) {
                    edges[i]->line_link = l1;
                    edges[i + 1]->line_link = l1;
                } else {
                    edges[i]->line_link = l2;
                    edges[i + 1]->line_link = l2;
                }

                if (i >= 4) {
                    auto l = std::make_shared<line>(*edges[i]->line_link);
                    l->a *= -1;
                    l->b *= -1;
                    l->c *= -1;
                    edges[i]->line_link = l;
                    edges[i + 1]->line_link = l;
                }
            }

            inf_node->e = edges[0];
            inner_vertex->e = edges[1];

            all_lines.push_back(l1);
            all_lines.push_back(l2);
        }

        DCEL(const std::vector<std::shared_ptr<line>> & lines)
            : max_vertex(0), max_edges(0)
        {
            auto edge_lines = std::vector<std::shared_ptr<line>>(3);
            auto vertices = std::vector<std::shared_ptr<vertex>>(3);
            auto edges = std::vector<std::shared_ptr<edge>>(6);

            line down(0, 1, 0), left(1, 0, 0), diag(1, 1, 0);
            edge_lines[0] = std::make_shared<line>(left);
            edge_lines[1] = std::make_shared<line>(down);
            edge_lines[2] = std::make_shared<line>(diag);

            std::vector<std::shared_ptr<line>> new_lines(lines);
            find_border_line(*edge_lines[0], 200, 1, new_lines);
            new_lines.push_back(edge_lines[0]);
            find_border_line(*edge_lines[1], 200, 1, new_lines);
            new_lines.push_back(edge_lines[1]);
            find_border_line(*edge_lines[2], -200, -1, new_lines);

            vertices[0] = std::make_shared<vertex>(edge_lines[0], edge_lines[2], max_vertex);
            vertices[1] = std::make_shared<vertex>(edge_lines[0], edge_lines[1], max_vertex);
            vertices[2] = std::make_shared<vertex>(edge_lines[1], edge_lines[2], max_vertex);

            for (int i = 0; i < 6; i++) {
                edges[i] = std::make_shared<edge>(max_edges);
                if (i % 2 == 0) {
                    edges[i]->hull_edge = true;
                }
            }

            for (int i = 0; i < 6; i += 2) {
                vertices[i / 2]->e = edges[(i + 4) % 6];

                edges[i]->origin = vertices[((i + 2) / 2) % 3];
                edges[i + 1]->origin = vertices[i / 2];

                edges[i]->twin = edges[i + 1];
                edges[i + 1]->twin = edges[i];

                edges[i]->line_link = edge_lines[i / 2];
                edges[i + 1]->line_link = edge_lines[i / 2];

                edges[i]->next = edges[(i + 4) % 6];
                edges[i + 1]->next = edges[(i + 3) % 6];

                edges[i]->prev = edges[(i + 2) % 6];
                edges[i + 1]->prev = edges[(i + 4) % 6];
            }

            inf_node = vertices[0];

            for (auto l : lines) {
                add_line_in_triangle(*l);
            }

            //add_line_in_triangle(*lines[0]);
        }

        void find_border_line(line & l, double d, int sign, const std::vector<std::shared_ptr<line>> & lines)
        {
            bool all_one_side = false;
            while (!all_one_side) {
                all_one_side = true;
                for (size_t i = 0; i < lines.size() - 1 && all_one_side; i++) {
                    for (size_t j = i + 1; j < lines.size() && all_one_side; j++) {
                        if (orientation_2d(lines[i]->a, lines[i]->b, lines[j]->a, lines[j]->b) == ZERO_DEAD) continue;
                        int line_sign = line_point_sign(l, *lines[i], *lines[j]);
                        all_one_side = line_sign * sign > 0;
                    }
                }
                if (!all_one_side) {
                    l.c += d;
                }
            }
        }

        void deep_copy_edge(int eid, const std::shared_ptr<edge> & other,
                            std::vector<std::shared_ptr<vertex>> & vertices,
                            std::vector<std::shared_ptr<edge>> & edges)
        {
            if (edges[eid]) return;

            edges[eid] = std::make_shared<edge>();
            edges[eid]->id = other->id;
            edges[eid]->line_link = other->line_link;
            edges[eid]->triangle_link = other->triangle_link;
            edges[eid]->hull_edge = other->hull_edge;
            edges[eid]->triangle_edge = other->triangle_edge;

            deep_copy_vertex(other->origin->id, other->origin, vertices, edges);
            deep_copy_edge(other->twin->id, other->twin, vertices, edges);
            deep_copy_edge(other->next->id, other->next, vertices, edges);
            deep_copy_edge(other->prev->id, other->prev, vertices, edges);

            edges[eid]->origin = vertices[other->origin->id];
            edges[eid]->twin = edges[other->twin->id];
            edges[eid]->next = edges[other->next->id];
            edges[eid]->prev = edges[other->prev->id];
        }

        void deep_copy_vertex(int vid, const std::shared_ptr<vertex> & other,
                              std::vector<std::shared_ptr<vertex>> & vertices,
                              std::vector<std::shared_ptr<edge>> & edges)
        {
            if (vertices[vid]) return;

            vertices[vid] = std::make_shared<vertex>();
            vertices[vid]->id = other->id;
            vertices[vid]->line1 = other->line1;
            vertices[vid]->line2 = other->line2;

            deep_copy_edge(other->e->id, other->e, vertices, edges);

            vertices[vid]->e = edges[other->e->id];
        }

        DCEL(const DCEL & other_dcel) // copy constructor, need for deleting vertices
        {
            max_vertex = other_dcel.max_vertex;
            max_edges = other_dcel.max_edges;
            std::vector<std::shared_ptr<vertex>> vertices(max_vertex);
            std::vector<std::shared_ptr<edge>> edges(max_edges);

            deep_copy_vertex(0, other_dcel.inf_node, vertices, edges);
            inf_node = vertices[0];

            all_lines.resize(other_dcel.all_lines.size());
            std::copy(other_dcel.all_lines.begin(), other_dcel.all_lines.end(), all_lines.begin());
        }

        bool edge_intersects_line(const line & l, const std::shared_ptr<edge> & e) const
        {
            if (is_ray(e)) {
                auto v = not_inf_vertex(e);
                return ray_line_intersection(l, *e->line_link, *v->line1, *v->line2);
            }

            auto v = e->origin;
            auto u = e->next->origin;
            return segment_line_intersection(l, *v->line1, *v->line2, *u->line1, *u->line2);
        }

        void add_line(const line & new_line)
        {
            std::shared_ptr<edge> inf_face_edge;
            auto e = inf_node->e;
            auto el = e->line_link;
            auto l = std::make_shared<line>(new_line);
            all_lines.push_back(l);

            if (orientation_2d(-el->b, el->a, -l->b, l->a) == NEG_DEAD) {
                inf_face_edge = e;
            } else {
                auto f = e->twin->next;
                auto fl = f->line_link;

                while (orientation_2d(-el->b, el->a, -l->b, l->a) == orientation_2d(-fl->b, fl->a, -l->b, l->a)) {
                    e = f;
                    f = f->twin->next;
                    el = e->line_link;
                    fl = f->line_link;
                }

                inf_face_edge = f;
            }

            auto crossed_edge = inf_face_edge;
            while (!edge_intersects_line(new_line, crossed_edge)) {
                crossed_edge = crossed_edge->next;
            }

            auto new_vertex = std::make_shared<vertex>(
                crossed_edge->line_link, l, max_vertex
            );

            auto line_edge1 = std::make_shared<edge>(max_edges); // line segment edges
            auto line_edge2 = std::make_shared<edge>(max_edges);
            auto part_edge1 = std::make_shared<edge>(max_edges); // partition half edges
            auto part_edge2 = std::make_shared<edge>(max_edges);

            new_vertex->e = part_edge1;

            part_edge1->origin = new_vertex;
            part_edge1->twin = crossed_edge->twin;
            part_edge1->next = (!crossed_edge->next->origin->line1 ? line_edge2 : crossed_edge->next);
            part_edge1->prev = line_edge2;
            part_edge1->line_link = crossed_edge->line_link;

            part_edge2->origin = new_vertex;
            part_edge2->twin = crossed_edge;
            part_edge2->next = crossed_edge->twin->next; // !!! may be next in line_edge1
            // part_edge2->prev -- cannot set now
            part_edge2->line_link = crossed_edge->line_link;

            line_edge1->origin = new_vertex;
            line_edge1->twin = line_edge2;
            line_edge1->next = inf_face_edge;
            line_edge1->prev = crossed_edge;
            line_edge1->line_link = l;

            line_edge2->origin = inf_node;
            line_edge2->twin = line_edge1;
            line_edge2->next = part_edge1;
            line_edge2->prev = (!crossed_edge->next->origin->line1 ? part_edge1 : inf_face_edge->prev);
            line_edge2->line_link = l;

            // redefine inf_node edge
            if (orientation_2d(-el->b, el->a, -l->b, l->a) == NEG_DEAD) {
                inf_node->e = line_edge2;
            }

            auto face_edge = crossed_edge->twin->next;
            if (crossed_edge->next->origin->line1) {
                crossed_edge->next->prev = part_edge1;
                inf_face_edge->prev->next = line_edge2;
            }
            inf_face_edge->prev = line_edge1;

            crossed_edge->next = line_edge1;
            crossed_edge->twin->twin = part_edge1;
            crossed_edge->twin = part_edge2;

            do {
                line_edge1 = std::make_shared<edge>(max_edges); // line segment edges
                line_edge2 = std::make_shared<edge>(max_edges);

                //line_edge1->origin = new_vertex;
                line_edge1->twin = line_edge2;
                line_edge1->next = part_edge2;
                //line_edge1->prev - cannot set now;
                line_edge1->line_link = l;

                line_edge2->origin = new_vertex;
                line_edge2->twin = line_edge1;
                //line_edge2->next - cannot set now;
                line_edge2->prev = part_edge1->twin;
                line_edge2->line_link = l;

                while (face_edge->id != part_edge1->twin->id && !edge_intersects_line(new_line, face_edge)) {
                    face_edge = face_edge->next;
                }

                if (face_edge->id == part_edge1->twin->id) { // end in last half plane
                    part_edge1->twin->next = part_edge2;

                    while (face_edge->origin->line1) {
                        face_edge = face_edge->next;
                    }
                    inf_face_edge = face_edge;

                    part_edge1->twin->next = line_edge2;
                    part_edge2->prev = line_edge1;

                    auto rev_line = std::make_shared<line>(new_line);
                    rev_line->a *= -1;
                    rev_line->b *= -1;
                    rev_line->c *= -1;

                    line_edge1->line_link = rev_line;
                    line_edge2->line_link = rev_line;

                    line_edge1->origin = inf_node;
                    line_edge1->prev = inf_face_edge->prev;
                    line_edge2->next = inf_face_edge;

                    inf_face_edge->prev->next = line_edge1;
                    inf_face_edge->prev = line_edge2;


                    break;
                }

                part_edge1->twin->next = line_edge2;
                part_edge2->prev = line_edge1;
                part_edge2->next->prev = part_edge2;

                auto new_vertex2 = std::make_shared<vertex>(
                    face_edge->line_link, l, max_vertex
                );

                auto new_part_edge1 = std::make_shared<edge>(max_edges); // partition half edges
                auto new_part_edge2 = std::make_shared<edge>(max_edges);

                line_edge1->prev = face_edge;
                line_edge2->next = new_part_edge1;
                line_edge1->origin = new_vertex2;

                new_vertex2->e = new_part_edge2;

                new_part_edge1->origin = new_vertex2;
                new_part_edge1->twin = face_edge->twin;
                new_part_edge1->next = face_edge->next;
                new_part_edge1->prev = line_edge2;
                new_part_edge1->line_link = face_edge->line_link;

                new_part_edge2->origin = new_vertex2;
                new_part_edge2->twin = face_edge;
                new_part_edge2->next = face_edge->twin->next; // remember to reset in last half plane
                //new_part_edge2->prev - cannot set now;
                new_part_edge2->line_link = face_edge->line_link;

                face_edge->next->prev = new_part_edge1;
                face_edge->twin->next->prev = new_part_edge2; // some problems here
                face_edge->next = line_edge1;
                //face_edge->twin->next = new_part_edge2;
                face_edge->twin->twin = new_part_edge1;
                face_edge->twin = new_part_edge2;

                crossed_edge = face_edge;
                face_edge = crossed_edge->twin->next;
                new_vertex = new_vertex2;
                part_edge1 = new_part_edge1;
                part_edge2 = new_part_edge2;
            } while (true);
        }

        void add_line_in_triangle(const line & new_line)
        {
            auto l = std::make_shared<line>(new_line);
            auto crossed_edge = inf_node->e;

            while (!edge_intersects_line(new_line, crossed_edge)) {
                crossed_edge = crossed_edge->next;
            }

            std::shared_ptr<vertex> new_vertex = std::make_shared<vertex>(
                crossed_edge->line_link, l, max_vertex
            );

            std::shared_ptr<edge> part_edge1 = std::make_shared<edge>(max_edges); // partition half edges
            std::shared_ptr<edge> part_edge2 = std::make_shared<edge>(max_edges);

            new_vertex->e = part_edge1;

            part_edge1->origin = new_vertex;
            part_edge1->twin = crossed_edge->twin;
            part_edge1->next = crossed_edge->next;
            part_edge1->prev = crossed_edge;
            part_edge1->line_link = crossed_edge->line_link;
            part_edge1->hull_edge = true;

            part_edge2->origin = new_vertex;
            part_edge2->twin = crossed_edge;
            part_edge2->next = crossed_edge->twin->next;
            //part_edge2->prev -- can't set now
            part_edge2->line_link = crossed_edge->line_link;

            auto face_edge = crossed_edge->twin->next;
            crossed_edge->next->prev = part_edge1;
            crossed_edge->next = part_edge1;
            crossed_edge->twin->next->prev = part_edge2;
            //crossed_edge->twin->next -- can't set now
            crossed_edge->twin->twin = part_edge1;
            crossed_edge->twin = part_edge2;

            do {
                std::shared_ptr<edge> line_edge1 = std::make_shared<edge>(max_edges); // line segment edges
                std::shared_ptr<edge> line_edge2 = std::make_shared<edge>(max_edges);

                part_edge1->twin->next = line_edge2;
                part_edge2->prev = line_edge1;

                // line_edge1->origin -- can't set now
                line_edge1->twin = line_edge2;
                line_edge1->next = part_edge2;
                //line_edge1->prev -- can't set now
                line_edge1->line_link = l;

                line_edge2->origin = new_vertex;
                line_edge2->twin = line_edge1;
                // line_edge2->next -- can't set noww
                line_edge2->prev = part_edge1->twin;
                line_edge2->line_link = l;

                while (face_edge->id != part_edge1->twin->id && !edge_intersects_line(new_line, face_edge)) {
                    face_edge = face_edge->next;
                }

                std::shared_ptr<vertex> new_vertex2 = std::make_shared<vertex>(
                    face_edge->line_link, l, max_vertex
                );

                std::shared_ptr<edge> new_part_edge1 = std::make_shared<edge>(max_edges); // partition half edges
                std::shared_ptr<edge> new_part_edge2 = std::make_shared<edge>(max_edges);

                line_edge1->origin = new_vertex2;
                line_edge1->prev = face_edge;
                line_edge2->next = new_part_edge1;

                new_vertex2->e = new_part_edge2;

                new_part_edge1->origin = new_vertex2;
                new_part_edge1->twin = face_edge->twin;
                new_part_edge1->next = face_edge->next;
                new_part_edge1->prev = line_edge2;
                new_part_edge1->line_link = face_edge->line_link;

                new_part_edge2->origin = new_vertex2;
                new_part_edge2->twin = face_edge;
                new_part_edge2->next = face_edge->twin->next;
                new_part_edge2->prev = new_part_edge1->twin; // this may change later
                new_part_edge2->line_link = face_edge->line_link;

                face_edge->next->prev = new_part_edge1;
                face_edge->twin->next->prev = new_part_edge2; // some problems here
                face_edge->next = line_edge1;
                face_edge->twin->next = new_part_edge2;
                face_edge->twin->twin = new_part_edge1;
                face_edge->twin = new_part_edge2;

                if (new_part_edge1->twin->hull_edge) { // end in last half plane
                    new_part_edge2->hull_edge = true;
                    break;
                }

                crossed_edge = face_edge;
                face_edge = crossed_edge->twin->next;
                new_vertex = new_vertex2;
                part_edge1 = new_part_edge1;
                part_edge2 = new_part_edge2;

            } while (true);
        }

        void get_intersected_edges(const line & new_line, std::vector<std::shared_ptr<edge>> & out) const
        {
            std::shared_ptr<edge> inf_face_edge;
            auto e = inf_node->e;
            auto el = e->line_link;
            auto l = std::make_shared<line>(new_line);

            if (orientation_2d(-el->b, el->a, -l->b, l->a) != POS_DEAD) {
                inf_face_edge = e;
            } else {
                auto f = e->twin->next;
                auto fl = f->line_link;

                while (true) {
                    dead_sign e_or = orientation_2d(-el->b, el->a, -l->b, l->a);
                    dead_sign f_or = orientation_2d(-fl->b, fl->a, -l->b, l->a);
                    if (e_or == ZERO_DEAD || (e_or != f_or && f_or != ZERO_DEAD)) break;
                    e = f;
                    f = f->twin->next;
                    el = e->line_link;
                    fl = f->line_link;
                }

                inf_face_edge = f;
            }

            auto crossed_edge = inf_face_edge;
            while (!edge_intersects_line(new_line, crossed_edge)) {
                crossed_edge = crossed_edge->next;
            }
            out.push_back(crossed_edge);

            auto face_edge = crossed_edge->twin->next;
            do {
                while (face_edge->id != crossed_edge->twin->id && !edge_intersects_line(new_line, face_edge)) {
                    face_edge = face_edge->next;
                }
                out.push_back(face_edge);

                if (face_edge->id == crossed_edge->twin->id) { // end in last half plane
                    break;
                }

                face_edge = face_edge->twin->next;
            } while (true);
        }

        orientation_t point2edge_orientation(const std::shared_ptr<edge> & e, const point_2 & c) const
        {
            if (is_ray(e)) {
                int res = line_position(*e->line_link, c);
                auto v = e->origin;
                if (res > 0) {
                    return (!v->line1
                            ? (is_direct_vector_right(*e->line_link) ? CG_RIGHT : CG_LEFT)
                            : (is_direct_vector_right(*e->line_link) ? CG_LEFT : CG_RIGHT));
                } else {
                    return (!v->line1
                            ? (is_direct_vector_right(*e->line_link) ? CG_LEFT : CG_RIGHT)
                            : (is_direct_vector_right(*e->line_link) ? CG_RIGHT : CG_LEFT));
                }
            }

            auto v = e->origin;
            auto u = e->next->origin;
            return point_segment_orientation(*v->line1, *v->line2, *u->line1, *u->line2, c);
        }


        std::shared_ptr<edge> get_face_by_point(const point_2 & p)
        {
            std::queue<std::shared_ptr<edge>> q;
            q.push(inf_node->e);
            std::vector<int> visited(max_edges);
            visited[inf_node->e->id] = 1;

            while (!q.empty()) {
                auto e = q.front();
                q.pop();

                if (visited[e->id] == 2) continue;
                visited[e->id] = 2;

                auto en = e;
                bool all_same = true;
                do {
                    visited[en->id] = 2;
                    all_same &= (CG_RIGHT != point2edge_orientation(en, p));
                    if (visited[en->twin->id] == 0) {
                        visited[en->twin->id] = 1;
                        q.push(en->twin);
                    }
                    en = en->next;
                } while (en->id != e->id);

                if (all_same) {
                    return e;
                }
            }

            return std::shared_ptr<edge>();
        }
    };
}
