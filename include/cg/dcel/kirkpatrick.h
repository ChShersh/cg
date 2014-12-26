#include <cg/dcel/dcel.h>

#include <set>

namespace cg {

    struct triangulation_level
    {
        DCEL dcel;
        std::vector< std::vector<std::shared_ptr<vertex>> > graph;
        std::vector< std::shared_ptr<vertex> > vertices;

        triangulation_level() {}

        triangulation_level(const DCEL & other_dcel, bool need_to_build_triangles) // building initial triangulation
            : dcel(other_dcel)
        {
            if (!need_to_build_triangles) return;

            std::queue<std::shared_ptr<vertex>> q;
            std::vector<bool> binded(dcel.max_edges);
            std::vector<int> visited(dcel.max_vertex);
            q.push(dcel.inf_node);
            visited[dcel.inf_node->id] = 1;

            while (!q.empty()) {
                std::shared_ptr<vertex> v = q.front();
                q.pop();

                std::shared_ptr<edge> e = v->e;
                do {
                    if (e->hull_edge || e->triangle_edge || binded[e->id]) {
                        binded[e->id] = true;
                        e = e->twin->next;
                        continue;
                    }
                    binded[e->id] = true;

                    auto last_edge = e;
                    auto f = e->next;
                    do {
                        binded[f->id] = true;

                        if (!visited[f->origin->id]) {
                            visited[f->origin->id] = 1;
                            q.push(f->origin);
                        }

                        std::shared_ptr<triangle_k> t = std::make_shared<triangle_k>(
                            line_cross(v->line1, v->line2),
                            line_cross(f->origin->line1, f->origin->line2),
                            line_cross(f->next->origin->line1, f->next->origin->line2)
                        );

                        f->triangle_link = std::make_shared<node_triangle>(t);
                        f->triangle_link->node_edge = f;
                        f->triangle_link->is_leaf = true;
                        f->triangle_link->depth = 0;

                        // join v and f->origin with triangle edge
                        auto next_f = f->next;
                        if (next_f->next->origin->id == v->id) {
                            f = next_f;
                            break;
                        }

                        std::shared_ptr<edge> tedge1 = std::make_shared<edge>(dcel.max_edges);
                        std::shared_ptr<edge> tedge2 = std::make_shared<edge>(dcel.max_edges);

                        tedge1->origin = v;
                        tedge1->twin = tedge2;
                        tedge1->next = next_f;
                        tedge1->prev = last_edge->prev;
                        tedge1->triangle_edge = true;

                        tedge2->origin = next_f->origin;
                        tedge2->twin = tedge1;
                        tedge2->next = last_edge;
                        tedge2->prev = f;
                        tedge2->triangle_edge = true;

                        last_edge->prev->next = tedge1;
                        last_edge->prev = tedge2;
                        f->next = tedge2;
                        next_f->prev = tedge1;

                        last_edge = tedge1;
                        f = next_f;
                    } while (f->next->origin->id != v->id);

                    binded[f->id] = true;

                    if (!visited[f->origin->id]) {
                        visited[f->origin->id] = 1;
                        q.push(f->origin);
                    }

                    e = e->twin->next;
                } while (e->id != v->e->id);
            }
        }

        void create_graph()
        {
            graph.clear();
            vertices.clear();
            graph.resize(dcel.max_vertex);
            vertices.resize(dcel.max_vertex);

            std::queue<std::shared_ptr<vertex>> que;
            std::vector<int> visited(dcel.max_vertex);
            visited[dcel.inf_node->id] = 1;
            que.push(dcel.inf_node);
            vertices[0] = dcel.inf_node;

            while (!que.empty()) {
                auto v = que.front();
                std::set<int> out_vertices;
                que.pop();

                std::shared_ptr<edge> e = v->e;
                do {
                    if (out_vertices.count(e->next->origin->id) == 0 &&
                        e->next->origin->id != v->id)
                    {
                        graph[v->id].push_back(e->next->origin);
                        out_vertices.insert(e->next->origin->id);
                    }

                    if (!visited[e->next->origin->id]) {
                        visited[e->next->origin->id] = 1;
                        que.push(e->next->origin);
                        vertices[e->next->origin->id] = e->next->origin;
                    }

                    e = e->twin->next;
                } while (e->id != v->e->id);
            }
        }
    };

    struct kirkpatrick_localization
    {
        DCEL dcel;
        DCEL hulled_dcel;
        std::vector<triangulation_level> levels;
        std::shared_ptr<node_triangle> root;
        int max_depth;

        kirkpatrick_localization() {}

        kirkpatrick_localization(const line & line1, const line & line2)
            : dcel(line1, line2)
        {}

        kirkpatrick_localization(const kirkpatrick_localization & other)
            : dcel(other.dcel)
        {}

        void add_line(const line & l)
        {
            dcel.add_line(l);
        }

        std::shared_ptr<edge> naive_localization(const point_2 & p)
        {
            return dcel.get_face_by_point(p);
        }

        std::shared_ptr<edge> fast_localization(const point_2 & p)
        {
            line l1(1, 0, -p.x), l2(0, 1, -p.y);
            if (!triangle_contains_convex_point(*root->t, l1, l2)) {
                return nullptr;
            }

            std::shared_ptr<node_triangle> node = root;
            while (!node->is_leaf) {
                for (auto child : node->children) {
                    if (triangle_contains_convex_point(*child->t, l1, l2)) {
                        node = child;
                        break;
                    }
                }
            }

            return node->node_edge;;
        }

        void build_triangulation(std::vector<std::vector<std::shared_ptr<vertex>>> & deleted_vertices)
        {
            deleted_vertices.clear();
            levels.clear();

            hulled_dcel = DCEL(dcel.all_lines);
            levels.push_back(triangulation_level(hulled_dcel, true));

            max_depth = 0;
            deleted_vertices.push_back(std::vector<std::shared_ptr<vertex>>());
            while (not_trivial_dcel(levels.back().dcel)) {
                max_depth++;
                levels.push_back(compress_level(levels.back(), deleted_vertices.back()));
                deleted_vertices.push_back(std::vector<std::shared_ptr<vertex>>());
            }
        }

        bool not_trivial_dcel(const DCEL & dcel)
        {
            std::queue<std::shared_ptr<vertex>> que;
            std::vector<int> visited(dcel.max_vertex);
            visited[dcel.inf_node->id] = 1;
            que.push(dcel.inf_node);
            int size = 1;

            while (!que.empty()) {
                auto v = que.front();
                que.pop();

                auto e = v->e;
                do {
                    if (!visited[e->next->origin->id]) {
                        visited[e->next->origin->id] = 1;
                        que.push(e->next->origin);
                        size++;
                    }

                    if (size > 3) return true;

                    e = e->twin->next;
                } while (e->id != v->e->id);
            }

            return size > 3;
        }

        triangulation_level compress_level(const triangulation_level & prev_level,
                                           std::vector<std::shared_ptr<vertex>> & deleted_vertices)
        {
            triangulation_level new_level(prev_level.dcel, false);
            new_level.create_graph();

            std::vector<int> marked(new_level.dcel.max_vertex);

            for (size_t i = 3; i < new_level.graph.size(); i++) { // ignore inf vertex, never delete it
                if (new_level.vertices[i] && new_level.graph[i].size() < 12 && !marked[i]) {
                    std::shared_ptr<vertex> del_v = new_level.vertices[i];
                    std::shared_ptr<edge> face_edge = nullptr;

                    std::vector<std::shared_ptr<node_triangle>> nodes;

                    deleted_vertices.push_back(del_v);

                    // delete vertex froc DCEL
                    auto e = del_v->e;
                    do {
                        if (e->hull_edge || e->twin->hull_edge) {
                            e = e->twin->next;
                            continue;
                        }

                        auto in_edge1 = e->twin->prev;
                        auto in_edge2 = e->next;
                        if (!face_edge) {
                            face_edge = in_edge2;
                        }

                        // get triangles of face
                        std::shared_ptr<node_triangle> t1 = find_triangle(in_edge1);
                        std::shared_ptr<node_triangle> t2 = find_triangle(in_edge2);
                        if (t1) nodes.push_back(t1);
                        if (t2) nodes.push_back(t2);

                        in_edge2->prev = in_edge1;
                        in_edge1->next = e->next;
                        if (in_edge2->origin->e->id == e->twin->id) {
                            in_edge2->origin->e = in_edge2;
                        }

                        auto v = in_edge1->origin;
                        auto u = in_edge2->origin;
                        auto t = in_edge2->next->origin;
                        if (precise_turn_predicate(*v->line1, *v->line2, *u->line1, *u->line2, *t->line1, *t->line2) == CG_COLLINEAR &&
                            vertex_size(e->next->origin) <= 2 && e->next->origin->id != 0
                            && in_edge1->id != in_edge2->twin->id)
                        {
                            if (face_edge->id == in_edge2->id) {
                                face_edge = in_edge1;
                            }
                            merge_two_edges(in_edge1, in_edge2);
                        }

                        e = e->twin->next;
                    } while (e->id != del_v->e->id);

                    if (e->hull_edge) {
                        del_v->e->twin->next = del_v->e->prev->twin;
                        del_v->e->prev->twin->prev = del_v->e->twin;
                    }

                    if (e->hull_edge && vertex_size(del_v) == 2) {
                        if (face_edge->id == e->prev->twin->id) {
                            face_edge = e->twin;
                        }
                        merge_two_edges(e->prev, e);
                    }

                    // retriangulate face
                    do {
                        auto v = face_edge->origin;
                        auto u = face_edge->next->origin;
                        auto s = face_edge->next->next->origin;

                        std::shared_ptr<triangle_k> triangle = std::make_shared<triangle_k>(
                            line_cross(v->line1, v->line2),
                            line_cross(u->line1, u->line2),
                            line_cross(s->line1, s->line2)
                        );
                        std::shared_ptr<node_triangle> node = std::make_shared<node_triangle>(triangle);
                        node->depth = max_depth;

                        if (is_triangle_face(face_edge)) {
                            face_edge->triangle_link = node;
                            node->node_edge = face_edge;
                            for (auto old_triangle : nodes) {
                                if (triangle_intersection(*triangle, *old_triangle->t)) {
                                    node->children.push_back(old_triangle);
                                }
                            }

                            if (!not_trivial_dcel(new_level.dcel)) {
                                root = node;
                            }

                            break;
                        }

                        bool is_ear = precise_turn_predicate(*v->line1, *v->line2, *u->line1, *u->line2, *s->line1, *s->line2) == CG_LEFT;
                        if (is_ear) {
                            auto e = face_edge;
                            do {
                                auto t = e->origin;
                                if (t->id != v->id && t->id != u->id && t->id != s->id) {
                                    is_ear = !triangle_contains_convex_point(*triangle, *t->line1, *t->line2);
                                }
                                e = e->next;
                            } while (is_ear && e->id != face_edge->id);
                        }

                        if (is_ear && !triangle_contains_star_point(*triangle, *del_v->line1, *del_v->line2)) {
                            // create new triangle edge
                            std::shared_ptr<edge> tedge1 = std::make_shared<edge>(new_level.dcel.max_edges);
                            std::shared_ptr<edge> tedge2 = std::make_shared<edge>(new_level.dcel.max_edges);

                            tedge1->origin = v;
                            tedge1->twin = tedge2;
                            tedge1->next = face_edge->next->next;
                            tedge1->prev = face_edge->prev;
                            tedge1->triangle_edge = true;

                            tedge2->origin = s;
                            tedge2->twin = tedge1;
                            tedge2->next = face_edge;
                            tedge2->prev = face_edge->next;
                            tedge2->triangle_edge = true;
                            tedge2->triangle_link = node;

                            face_edge->prev->next = tedge1;
                            face_edge->prev = tedge2;
                            face_edge->next->next->prev = tedge1;
                            face_edge->next->next = tedge2;

                            for (auto old_triangle : nodes) {
                                if (triangle_intersection(*triangle, *old_triangle->t)) {
                                    node->children.push_back(old_triangle);
                                }
                            }

                            node->node_edge = face_edge;
                            face_edge = tedge1;
                        } else {
                            face_edge = face_edge->next;
                        }
                    } while (true);

                    // mark neighbour nodes as visited
                    for (auto u : new_level.graph[i]) {
                        marked[u->id] = 1;
                    }
                }
                marked[i] = 1;
            }

            return new_level;
        }

        int vertex_size(const std::shared_ptr<vertex> & v)
        {
            int size = 0;
            auto e = v->e;
            do {
                size++;
                e = e->twin->next;
            } while (e->id != v->e->id);
            return size;
        }

        bool is_triangle_face(const std::shared_ptr<edge> & e)
        {
            auto f = e;
            int size = 0;
            do {
                size++;
                if (size > 3) return false;
                f = f->next;
            } while (f->id != e->id);
            return size <= 3;
        }

        void merge_two_edges(std::shared_ptr<edge> in_edge1, std::shared_ptr<edge> in_edge2)
        {
            in_edge2->next->prev = in_edge1;
            in_edge1->next = in_edge2->next;
            in_edge1->twin->next->prev = in_edge2->twin;
            in_edge2->twin->next = in_edge1->twin->next;
            in_edge1->twin = in_edge2->twin;
            in_edge2->twin->twin = in_edge1;
        }

        std::shared_ptr<node_triangle> find_triangle(const std::shared_ptr<edge> & e)
        {
            auto f = e;
            for (int i = 0; i < 3; i++) {
                if (!f->triangle_link) {
                    f = f->next;
                } else {
                    auto t = f->triangle_link;
                    f->triangle_link = std::shared_ptr<node_triangle>();
                    return t;
                }
            }

            return std::shared_ptr<node_triangle>();
        }
    };
}
