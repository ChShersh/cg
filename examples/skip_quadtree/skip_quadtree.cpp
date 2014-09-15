#include <QColor>
#include <QApplication>

#include <cg/visualization/viewer_adapter.h>
#include <cg/visualization/draw_util.h>

#include <cg/io/point.h>

#include <cg/trees/quadtree.h>
#include <cg/trees/compressed_quadtree.h>
#include <cg/trees/skip_quadtree.h>

#include <misc/random_utils.h>

#include <set>

using cg::point_2f;
using cg::point_2;
using cg::rectangle_2;
using cg::range;

using cg::quadtree;
using cg::compressed_quadtree;
using cg::QuadNode;
using cg::skip_quadtree;

struct quad_tree_viewer : cg::visualization::viewer_adapter
{
   quad_tree_viewer()
      : quadtree_type(quad), mode(inserting), lx(-225), ly(-225), rx(225), ry(225),
        current_layer(0), is_found(false), corners(0), eps(10)
   {
       quad_tree = new quadtree<double>(lx, ly, rx, ry);
       compressed_tree = compressed_quadtree<double>(lx, ly, rx, ry);
       skip_tree = skip_quadtree<double>(lx, ly, rx, ry);
   }

   void draw_rectangle(const point_2 & dlc, const point_2 & urc,
                       cg::visualization::drawer_type & drawer, QColor color) const
   {
       drawer.set_color(color);
       drawer.draw_line(point_2{dlc.x, dlc.y}, point_2{dlc.x, urc.y});
       drawer.draw_line(point_2{dlc.x, urc.y}, point_2{urc.x, urc.y});
       drawer.draw_line(point_2{urc.x, urc.y}, point_2{urc.x, dlc.y});
       drawer.draw_line(point_2{urc.x, dlc.y}, point_2{dlc.x, dlc.y});
   }

   void draw_rectangle(const rectangle_2 & rect, cg::visualization::drawer_type & drawer, QColor color) const
   {
       draw_rectangle(point_2{rect.x.inf, rect.y.inf}, point_2{rect.x.sup, rect.y.sup}, drawer, color);
   }

   void visualize_quadtree(quadtree<double> * tr, cg::visualization::drawer_type & drawer) const
   {
      if (tr->is_leaf) {
         if (tr->point != boost::none) {
            drawer.set_color(Qt::white);
            drawer.draw_point(tr->point.get(), 2);
         }
      } else {
         drawer.set_color(Qt::green);
         drawer.draw_line(point_2{tr->lx, (tr->ly + tr->ry) / 2}, point_2{tr->rx, (tr->ly + tr->ry) / 2});
         drawer.draw_line(point_2{(tr->lx + tr->rx) / 2, tr->ly}, point_2{(tr->lx + tr->rx) / 2, tr->ry});

         for (int i = 0; i < 4; i++) {
             visualize_quadtree(tr->children[i], drawer);
         }
      }
   }

   void visualize_compressed_quadtree(const std::shared_ptr<QuadNode<double>> & tr,
                                      cg::visualization::drawer_type & drawer) const
   {
       if (tr->is_leaf) {
           if (tr->point != boost::none) {
               drawer.set_color(Qt::white);
               drawer.draw_point(tr->point.get(), 2);
           }
       } else {
           int non_empty = 0;
           for (int i = 0; i < 4; i++) {
               if (tr->children[i]) {
                   auto ch = tr->children[i];
                   non_empty++;

                   draw_rectangle(point_2{ch->lx, ch->ly}, point_2{ch->rx, ch->ry}, drawer, Qt::green);
                   visualize_compressed_quadtree(ch, drawer);
               }
           }

//           if (non_empty > 1) {
//               for (int i = 0; i < 4; i++) {
//                   if (tr->children[i]) {
//                       auto c = tr->coordinates_by_id(i);
//                       auto ax = c[0], ay = c[1], bx = c[2], by = c[3];
//
//                       draw_rectangle(point_2{ax, ay}, point_2{bx, by}, drawer, Qt::green);
//                   }
//               }
//           }
       }
   }

   void draw(cg::visualization::drawer_type & drawer) const
   {

      draw_rectangle(point_2{lx, ly}, point_2{rx, ry}, drawer, Qt::green);

      drawer.set_color(Qt::white);
      for (auto p : outside) {
          drawer.draw_point(p, 2);
      }

      switch (quadtree_type) {
          case quad : visualize_quadtree(quad_tree, drawer); break;
          case comp : visualize_compressed_quadtree(compressed_tree.root, drawer); break;
          case skip : visualize_compressed_quadtree(skip_tree.trees[current_layer].root, drawer); break;
      }

      if (is_found) {
          switch (quadtree_type) {
              case quad :
                  draw_rectangle(point_2{fq->lx, fq->ly}, point_2{fq->rx, fq->ry}, drawer, Qt::blue); break;
              case comp :
                  draw_rectangle(point_2{fc->lx, fc->ly}, point_2{fc->rx, fc->ry}, drawer, Qt::blue); break;
              case skip :
              {
                  auto fn = fs[current_layer];
                  draw_rectangle(point_2{fn->lx, fn->ly}, point_2{fn->rx, fn->ry}, drawer, Qt::blue);
                  break;
               }
          }
      }

      if (corners > 0) {
          auto dlc = point_2{std::min(rlc.x, rrc.x), std::min(rlc.y, rrc.y)};
          auto urc = point_2{std::max(rlc.x, rrc.x), std::max(rlc.y, rrc.y)};

          auto rect = rectangle_2(dlc, urc);
          auto eps_rect = rectangle_2{
              range{rect.x.inf - eps, rect.x.sup + eps},
              range{rect.y.inf - eps, rect.y.sup + eps}
          };

          draw_rectangle(rect, drawer, Qt::darkMagenta);
          draw_rectangle(eps_rect, drawer, Qt::cyan);

          std::vector<point_2> ans;
          switch (quadtree_type) {
              case quad : quad_tree->rectangle_query(rect, eps, ans); break;
              case comp : compressed_tree.rectangle_query(rect, eps, ans); break;
              case skip : skip_tree.approx_rect_query(rect, eps, ans, current_layer); break;
          }

          drawer.set_color(Qt::blue);
          for (auto pt : ans) {
              drawer.draw_point(pt, 5);
          }
      }

      if (cur_cursor != boost::none) {
          drawer.set_color(Qt::yellow);
          drawer.draw_point(cur_cursor.get(), 5);
      }
   }

   void print(cg::visualization::printer_type & p) const
   {
      p.corner_stream() << "press mouse rbutton to add point, "
                        << "press double-click to clear" << cg::visualization::endl
                        << "press 'i' for insert mode, 'd' for delete mode, 'f' for find mode, 'r' for rectangle mode" << cg::visualization::endl
                        << "press 'q', 'c' or 's' to change tree structure; press 'n' to add random point from square" << cg::visualization::endl
                        << "skip quadtree: bottom_left = " << point_2{lx, ly} << ", "
                        <<                "top_right = "   << point_2{rx, ry} << cg::visualization::endl << cg::visualization::endl
                        << "current mode : ";

      switch (mode)
      {
          case inserting      : p.corner_stream() << "Insert" << cg::visualization::endl; break;
          case deleting       : p.corner_stream() << "Delete" << cg::visualization::endl; break;
          case finding        : p.corner_stream() << "Find" << cg::visualization::endl; break;
          case rectangle_mode : p.corner_stream() << "Rectangle" << cg::visualization::endl; break;
      }

      switch (quadtree_type)
      {
         case quad : p.corner_stream() << "using: quadtree" << cg::visualization::endl; break;
         case comp : p.corner_stream() << "using: compressed quadtree" << cg::visualization::endl; break;
         case skip : p.corner_stream() << "using: skip-quadtree; press left-right arrows to change layer: "
                                       << current_layer + 1 << "(" << skip_tree.trees.size() << ")" << cg::visualization::endl; break;
      }
   }

   void add_point(const point_2 & p)
   {
       if (p.x < lx || rx <= p.x || p.y < ly || ry <= p.y) {
           outside.insert(p);
           return;
       }

       switch (quadtree_type)
       {
           case quad :
               quad_tree->insert(p);
               quadtree_points.insert(p);
               break;
           case comp :
               compressed_tree.insert(p);
               compressed_points.insert(p);
               break;
           case skip :
               skip_tree.insert(p);
               skip_points.insert(p);
               is_found = false;
               break;
       }
   }

   void remove_point(const point_2 & pt)
   {
       switch (quadtree_type)
       {
           case quad : set_cur_cursor(pt, quadtree_points); break;
           case comp : set_cur_cursor(pt, compressed_points); break;
           case skip : set_cur_cursor(pt, skip_points); break;
       }

       if (cur_cursor == boost::none) {
           return;
       }

       auto p = cur_cursor.get();

       if (p.x < lx || rx <= p.x || p.y < ly || ry <= p.y) {
           outside.erase(p);
       } else {
           switch (quadtree_type)
           {
               case quad :
                   quad_tree->remove(p);
                   quadtree_points.erase(p);
                   break;
               case comp :
                   compressed_points.erase(p);
                   break;
               case skip :
                   skip_points.erase(p);
                   is_found = false;
                   break;
           }
       }
   }

   void find_point(const point_2 & p)
   {
       is_found = !(p.x < lx || rx <= p.x || p.y < ly || ry <= p.y);
       switch (quadtree_type) {
           case quad : fq = quad_tree->find(p); break;
           case comp : fc = compressed_tree.find(p); break;
           case skip : fs = skip_tree.search_all_levels(p); break;
       }
   }

   bool on_release(const point_2f & p)
   {
      switch (mode) {
          case inserting      : add_point(p); break;
          case deleting       : remove_point(p); break;
          case finding        : find_point(p); break;
          case rectangle_mode :
              if (corners == 0) {
                  rlc = p;
                  rrc = p;
                  corners = 1;
              } else if (corners == 1) {
                  rrc = p;
                  corners = 2;
              } else {
                  corners = 0;
              }
              break;
      }

      return true;
   }

   bool on_double_click(const point_2f & p)
   {
       switch (quadtree_type)
       {
           case quad :
              delete quad_tree;
              quad_tree = new quadtree<double>(lx, ly, rx, ry);
              quadtree_points.clear();
              break;
           case comp :
              compressed_tree = compressed_quadtree<double>(lx, ly, rx, ry);
              compressed_points.clear();
              break;
           case skip :
              skip_tree = skip_quadtree<double>(lx, ly, rx, ry);
              skip_points.clear();
              break;
       }
       outside.clear();
       current_layer = 0;
       corners = 0;
       is_found = false;

       return true;
   }

   bool on_key(int key)
   {
       switch (key)
       {
           case Qt::Key_I     : mode = inserting; is_found = false; break;
           case Qt::Key_D     : mode = deleting; is_found = false; break;
           case Qt::Key_F     : mode = finding; break;
           case Qt::Key_R     : if (corners > 0) corners = 2; mode = rectangle_mode; is_found = false; break;
           case Qt::Key_Q     : quadtree_type = quad; is_found = false; break;
           case Qt::Key_C     : quadtree_type = comp; is_found = false; break;
           case Qt::Key_S     : quadtree_type = skip; is_found = false; break;
           case Qt::Key_Left  :
               if (quadtree_type == skip)
                   current_layer = std::max(current_layer - 1, 0);
               break;
           case Qt::Key_Right :
               if (quadtree_type == skip)
                   current_layer = std::min(current_layer + 1, (int) (skip_tree.trees.size() - 1));
               break;
           case Qt::Key_N     :
           {
               using util::random_point;
               auto p = random_point(lx, rx, ly, ry);
               add_point(p);
               break;
           }
           default            : return false;
       }

       return true;
   }

   bool set_cur_cursor(const point_2 & p, const std::set<point_2> & all_points)
   {
       cur_cursor = boost::none;
       double max_r = 0;

       for (auto it = all_points.begin(); it != all_points.end(); it++) {
           double current_r = (p.x - it->x) * (p.x - it->x) + (p.y - it->y) * (p.y - it->y);

           if ((cur_cursor == boost::none && current_r < 100) || (cur_cursor != boost::none && current_r < max_r)) {
               cur_cursor = *it;
               max_r = current_r;
           }
       }

       if (cur_cursor == boost::none) {
           for (auto it = outside.begin(); it != outside.end(); it++) {
               double current_r = (p.x - it->x) * (p.x - it->x) + (p.y - it->y) * (p.y - it->y);

               if ((cur_cursor == boost::none && current_r < 100) || (cur_cursor != boost::none && current_r < max_r)) {
                   cur_cursor = *it;
                   max_r = current_r;
               }
           }
       }

       return cur_cursor != boost::none;
   }

   bool on_move(const point_2f & p)
   {
       switch (quadtree_type)
       {
           case quad : set_cur_cursor(p, quadtree_points); break;
           case comp : set_cur_cursor(p, compressed_points); break;
           case skip : set_cur_cursor(p, skip_points); break;
       }

       if (mode == rectangle_mode && corners == 1) {
           rrc = p;
       }
       return true;
   }

private:
   enum tree_type
   {
       quad = 1,
       comp = 2,
       skip = 3
   };

   enum cursor_mode
   {
       inserting = 0,
       deleting = 1,
       finding = 2,
       rectangle_mode = 3
   };

   // trees
   tree_type quadtree_type;

   quadtree<double> * quad_tree;
   std::set<point_2> quadtree_points;

   mutable compressed_quadtree<double> compressed_tree;
   std::set<point_2> compressed_points;

   mutable skip_quadtree<double> skip_tree;
   int current_layer;
   std::set<point_2> skip_points;

   cursor_mode mode;
   boost::optional<point_2> cur_cursor;
   std::set<point_2> outside;

   bool is_found;
   quadtree<double> * fq;
   std::shared_ptr<QuadNode<double>> fc;
   std::vector<std::shared_ptr<QuadNode<double>>> fs;

   double lx, ly, rx, ry, eps;
   point_2 rlc, rrc;
   int corners = 0;
};

int main(int argc, char ** argv)
{
   QApplication app(argc, argv);
   quad_tree_viewer viewer;
   cg::visualization::run_viewer(&viewer, "quadtrees");
}
