#pragma once

const int max_num_results = 100;

template <typename Tv>
__inline__ __host__ __device__ Tv _square_dist(Tv qx, Tv qy, Tv qz, Tv gx, Tv gy,
                                        Tv gz) {
                return pow(qx - gx, 2) + pow(qy - gy, 2) + pow(qz - gz, 2);
}

template <template <typename> class Array, typename Tvo, typename Tvi>
class ArgNearest {

        size_t num_points, point_dim;
        
        Tvo *out;
        size_t out_line;

        const Tvi *query;
        size_t query_line;
        
        const Tvi *points;
        size_t points_line;

        public:
         ArgNearest(Array<Tvo>& out_, const Array<Tvi>& query_,
                    const Array<Tvi>& points_) {
                 num_points = points_.nx;
                 point_dim = points_.ny;

                 out = out_();
                 out_line = out_.line;

                 query = query_();
                 query_line = query_.line;

                 points = points_();
                 points_line = points_.line;
         }

         VOID_OPERATOR() const {
                 int idx[max_num_results];

                 for (int k = 0; k <= yj; ++k) {
                         idx[k] = -1;
                        Tvi min_dist = -1;
                 for (int i = 0; i < num_points; ++i) {
                         // Skip previous minima
                         bool is_seen = false;
                         for (int j = 0; j <= k; ++j) {
                                 if (i == idx[j]) {
                                         is_seen = true;
                                         break;
                                 }
                         }
                         if (is_seen) continue;

                         Tvi qx = 0.0;
                         Tvi qy = 0.0;
                         Tvi qz = 0.0;
                         Tvi px = 0.0;
                         Tvi py = 0.0;
                         Tvi pz = 0.0;

                         if (point_dim >= 1) { 
                                qx = query[xi];
                                px = points[i];
                         }
                         if (point_dim >= 2) { 
                                qy = query[xi + query_line];
                                py = points[i + points_line];
                         }
                         if (point_dim == 3) {
                                qz = query[xi + 2 * query_line];
                                pz = points[i + 2 * points_line];
                         }

                         Tvi dist = _square_dist(qx, qy, qz, px, py, pz);

                         if (min_dist == -1 || dist < min_dist) {
                                 idx[k] = i;
                                 min_dist = dist;
                         }
                 }
                 }

                 size_t pos = xi + out_line * yj;
                 out[pos] = (size_t)idx[yj];
        }

};
