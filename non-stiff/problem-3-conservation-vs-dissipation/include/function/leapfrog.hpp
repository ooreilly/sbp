#pragma once

template <typename pFunctor, typename vFunctor, typename callbackFunctor, typename Tv=double>
void leapfrog(pFunctor& pf, vFunctor& vf, callbackFunctor& callback,
              const Tv T, const Tv dt) {

        callback(0, 0, dt);
        int num_steps = (int) ((T - 0.5 * dt ) / dt);
        for (int step = 0; step < num_steps; ++step) {
                Tv t = step * dt;
                pf(t, dt);
                vf(t + 0.5 * dt, dt);
                callback(step + 1, t + dt, dt);
        }

        // Take final correction step to reach the desired final time
        Tv new_dt = T - (num_steps) * dt;
        Tv t = num_steps * dt;
        if (new_dt > 0) {
                pf(t, dt);
                vf(t + 0.5 * dt, dt);
        }
        callback(num_steps, T, dt, 1);
}
