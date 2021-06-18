#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <unistd.h>
#include <pwd.h>
#include <cmath>
#include <boost/property_tree/ptree.hpp>

using namespace std;


bool intersect(double ax1, double ay1, double ax2, double ay2, double bx1, double by1, double bx2, double by2) {
    // segment_a = [[ax1, ay1], [ax2, ay2]], segment_b = [[bx1, by1], [bx2, by2]]
    double Aa, Ab, Ba, Bb, x_intersec;
    // Check for overlap: If the largest x coordinate of segment a is smaller than the smallest x coordinate
    // of segment b then there can be no intersection. (same for y)
    if (max(ax1, ax2) < min(bx1, bx2)){
        return false;
    }
    if (max(ay1, ay2) < min(by1, by2)){
        return false;
    }

    // If the is a mutual interval calculate the x coordinate of that intersection point and check if it is in the interval.
    // Calculate fa(a) = Aa*x + Ba = y and fb(x) = Ab*x + Bb = y
    Aa = (ay1 - ay2) / (ax1 - ax2);  // slope of segment a
    Ab = (by1 - by2) / (bx1 - bx2);  // slope of segment b
    Ba = ay1 - Aa * ax1;  // y intercept of segment a
    Bb = by1 - Ab * bx1;  // y intercep of segment b
    x_intersec = (Bb - Ba) / (Aa - Ab);  // x coordinate of intersection point
    if (x_intersec < max(min(ax1, ax2), min(bx1, bx2)) or x_intersec > min(max(ax1, ax2), max(bx1, bx2))){
        return false;
    }
    else{
        return true;
    }
}


class Adaptive_leaky_integrate_and_fire {
private:
    double _v, _a;
    float _mu, _d, _taua, _deltaa;
    double _isi;

    double _xi = 0;
    double _dt = 0.0001;

    std::random_device rd;
    std::mt19937 _generator;
    std::normal_distribution<double> dist;

public:
    Adaptive_leaky_integrate_and_fire(double v0, double a0, float mu, float D, float tau_a, float delta_a)
            : _generator(rd()),
              dist(std::normal_distribution<double>(
                      0, 1)) {
        _v = v0;
        _a = a0;
        _mu = mu;
        _d = D;
        _taua = tau_a;
        _deltaa = delta_a;
        _isi = period();

    }

    bool integrate_stochastic(){
        _xi = dist(_generator);
        _v += (_mu  - _v - _a) * _dt + sqrt(2 * _d * _dt) * _xi;
        _a += (-_a/_taua) * _dt;
        if (_v > 1.0){
            _v = 0;
            _a += _deltaa;
            return true;
        }
        else{
            return false;
        }
    }

    bool integrate_deterministic(){
        _v += (_mu  - _v - _a) * _dt;
        _a += (-_a/_taua) * _dt;
        if (_v > 1.0){
            _v -= 1;
            _a += _deltaa;
            return true;
        }
        else{
            return false;
        }
    }

    vector<double> integrate_deterministic(double v, double a) const{
        v += (_mu  - v - a) * _dt;
        a += (-a/_taua) * _dt;
        if (v > 1.0){
            v -= 1;
            a += _deltaa;
        }
        return {v, a};

    }

    vector<vector<double>> get_limit_cycle(){
        int spike_count;
        bool fired;

        spike_count = 0;
        // Allow the adaptive LIF model to evolve for some time so that we are close to the limit cycle
        while(spike_count < 100){
            fired = integrate_deterministic();
            if (fired){
                spike_count += 1;
            }
        }
        // Integrate one more time (from reset to threshold) and save as limit cycle
        std::vector<std::vector<double>> limit_cycle;
        std::vector<double> p;
        while(true){
            p = {_v, _a};
            limit_cycle.push_back(p);
            fired = integrate_deterministic();
            if (fired){
                return limit_cycle;
            }
        }
    }

    double period(){
        double v, a, t;
        int count_spikes;
        bool fired;

        _v=0;
        _a=0;
        count_spikes = 0;
        while(count_spikes < 100){
             fired = integrate_deterministic();
             if(fired){
                 count_spikes += 1;
             }
        }
        t = 0;
        count_spikes = 0;
        while (count_spikes < 100){
            fired = integrate_deterministic();
            t += _dt;
            if(fired){
                count_spikes += 1;
            }
        }
        return t /count_spikes;
    }

    vector<double> integrate_backwards_for_T(double v0, double a0) const{
        double v = v0;
        double a = a0;
        double t = 0;
        while(t < _isi) {
            t += _dt;
            v -= (_mu - v - a) * _dt;
            a -= (-a / _taua) * _dt;
        }
        return {v, a};
    }

    double get_deltaa() const{
        return _deltaa;
    }

    double get_dt() const{
        return _dt;
    }

    double get_a() const{
        return _a;
    }

    double get_v() const{
        return _v;
    }

    void set_a(double a){
        _a = a;
    }

    void set_v(double v){
        _v = v;
    }

};


class Bob_the_isochrone_builder{
private:
    Adaptive_leaky_integrate_and_fire& _model;
    double _isi;
    double _v0, _a0;
    double _dt, _deltaa;
    double _step_width;
    vector<double> _start_point_lower_branch, _start_point_upper_branch;
    vector<vector<double>> _isochrone;

public:
    Bob_the_isochrone_builder(Adaptive_leaky_integrate_and_fire& model, double isi, double v0, double a0, double step_width) : _model(model){
        _isi = isi;
        _v0 = v0;
        _a0 = a0;
        _dt = _model.get_dt();
        _deltaa = _model.get_deltaa();
        _step_width = step_width;
        _isochrone.push_back({v0, a0});
    }

    void insert_to_isochrone(int idx, double v0, double a0){
        _isochrone.insert(_isochrone.begin() + idx, {v0, a0});
    }

    vector<double> get_start_point_lower_branch() const{
        return _start_point_lower_branch;
    }

    vector<double> get_start_point_upper_branch() const{
        return _start_point_upper_branch;
    }

    vector<vector<double>> get_isochrone(){
        return _isochrone;
    }

    vector<double> next_point_right(vector<vector<double>>& list) const{
        vector<double> p1, p0;
        double v1, a1, v0, a0, dv, da, v_cand, a_cand, dif;
        int n =  list.size();
        if(list.size() > 1){
            p1 = list[n-1];
            p0 = list[n-2];
            v1 = p1[0];
            a1 = p1[1];
            v0 = p0[0];
            a0 = p0[1];
            dv = v1 - v0;
            da = a1 - a0;
            da = _step_width * (da / abs(dv));
            dv = _step_width * (dv / abs(dv));
            v_cand = v1 + dv;
            a_cand = a1 + da;
        }
        else{
            p1 = list[n-1];
            v1 = p1[0];
            a1 = p1[1];
            dif = fmod(v1, 0.05);
            v_cand = v1 + 0.05 - dif;
            a_cand = a1;
        }
        return {v_cand, a_cand};
    }

    vector<double> next_point_left(vector<vector<double>>& list) const{
        vector<double> p1, p0;
        double v1, a1, v0, a0, dv, da, v_cand, a_cand, dif;
        if(list.size() > 1){
            p1 = list[0];
            p0 = list[1];
            v1 = p1[0];
            a1 = p1[1];
            v0 = p0[0];
            a0 = p0[1];
            dv = v1 - v0;
            da = a1 - a0;
            da = _step_width * (da / abs(dv));
            dv = _step_width * (dv / abs(dv));
            v_cand = v1 + dv;
            a_cand = a1 + da;
        }
        else{
            p1 = list[0];
            v1 = p1[0];
            a1 = p1[1];
            dif = fmod(v1, 0.05);
            if(dif <= 0.0001){
                dif = 0.05;
            }
            v_cand = v1 - dif;
            a_cand = a1;
        }
        return {v_cand, a_cand};
    }

    double return_time_point_isochrone(double v, double a) {
        double t;
        vector<double> p, p_before, p_after;
        double v_before, a_before, v_after, a_after;
        double v_iso_before, a_iso_before, v_iso_after, a_iso_after;
        bool passed_isochrone, fired;
        int count_ref;

        t = 0;
        count_ref = 0;
        // Let the point (v, a) evolve until it hits the presumed isochrone
        while (true) {
            v_before = v;
            a_before = a;
            p = _model.integrate_deterministic(v, a);
            v = p[0];
            a = p[1];
            if(abs(v - v_before) > 0.5){
                v_after = v + 1.;
                a_after = a - _deltaa;
            }
            else {
                v_after = v;
                a_after = a;
            }
            t += _dt;
            count_ref += 1;
            // For every segment of the isochrone check if this segment and the segment that describes
            //  the change of v, a ((v_tmp, a_tmp), (v, a)) intersect.
            for(int i=0; i < _isochrone.size()-1; i++) {
                p_before = _isochrone[i];
                p_after = _isochrone[i + 1];
                v_iso_before = p_before[0];
                a_iso_before = p_before[1];
                v_iso_after = p_after[0];
                a_iso_after = p_after[1];
                if (std::abs(v_iso_after - v_iso_before) > 0.5) {
                    continue;
                }
                passed_isochrone = intersect(v_before, a_before, v_after, a_after, v_iso_before, a_iso_before,
                                             v_iso_after, a_iso_after);
                if (passed_isochrone and count_ref > 1) {
                    return t;
                }
            }
        }
    }

    vector<double> check_point_on_isochrone(double v_candidate, double a_candidate){
        cout << "Start new branch" << endl;
        bool on_isochrone;
        double T;

        _isochrone.push_back({v_candidate, a_candidate});
        on_isochrone = false;
        while(!on_isochrone){
            T = return_time_point_isochrone(v_candidate, a_candidate);
            if( abs(T - _isi) / _isi > 0.002){
                if(T > _isi){
                    cout << std::setprecision(4) << "ISI too long " <<  abs(T - _isi) / _isi << " " << "v: " << v_candidate << " a: "
                         << a_candidate << endl;
                    a_candidate = a_candidate - max(0.001, 0.1 * abs(T - _isi) / _isi);
                }
                else{
                    cout << std::setprecision(4) << "ISI too short " << abs(T - _isi) / _isi << " " << "v: " << v_candidate << " a: "
                         << a_candidate << endl;
                    a_candidate = a_candidate + max(0.001, 0.1 * abs(T - _isi) / _isi);
                }
                _isochrone.back() = {v_candidate, a_candidate};
            }
            else{
                cout << "on point" << endl;
                on_isochrone = true;
            }
        }
        return {v_candidate, a_candidate};
    }

    vector<vector<double>> construct_isochrone_branch_left_from_point(double v0, double a0, double vres, double vmin, int insert_point) {
        // Add initial points (v0, a0) to left part of isochrone
        vector<double> p;
        double v_candidate, a_candidate, T;
        bool on_isochrone;
        vector<vector<double>> isochrone_left = {{v0, a0}};

        cout << "Start left branch of isochrone" << endl;
        while (true) {
            p = next_point_left(isochrone_left);
            v_candidate = p[0];
            a_candidate = p[1];
            double a_before_on_isochrone = isochrone_left[0][1];

            if (v_candidate > vmin and a_before_on_isochrone > 0) {
                // Insert (v0, a0) to the general isochrone this is important because i calculate the return time to this
                // isochrone. Note that (v0, a0) is added at the beginning of isochrone.
                _isochrone.insert(_isochrone.begin() + insert_point, {v_candidate, a_candidate});
                on_isochrone = false;
                // Calculate return time from (v0, a0) to isochrone and adjust (v0, a0) until |T - ISI| / ISI < 0.0005
                while (!on_isochrone) {
                    T = return_time_point_isochrone(v_candidate, a_candidate);
                    if (abs(T - _isi) / _isi > 0.001) {
                        if (T > _isi) {
                            cout << std::setprecision(4) << "ISI too long " << abs(T - _isi) / _isi << " " << "v: " << v_candidate << " a: "
                                 << a_candidate << endl;
                            a_candidate = a_candidate - max(0.001, 0.1 * abs(T - _isi) / _isi);
                        } else {
                            cout << std::setprecision(4) << "ISI too short " << abs(T - _isi) / _isi << " " << "v: " << v_candidate << " a: "
                                 << a_candidate << endl;
                            a_candidate = a_candidate + max(0.001, 0.1 * abs(T - _isi) / _isi);
                        }
                        _isochrone[insert_point] = {v_candidate, a_candidate};
                    }
                    else {
                        cout << "on point" << endl;
                        on_isochrone = true;
                        if (abs(v_candidate - vres) < 0.001) {
                            _start_point_lower_branch = {1, a_candidate - _deltaa};
                        }
                        isochrone_left.insert(isochrone_left.begin() + 0, {v_candidate, a_candidate});
                    }
                }
            }
            else {
                return isochrone_left;
            }
        }
    }

    vector<vector<double>> construct_isochrone_branch_right_from_point(double v0, double a0, double vmax){
        vector<double> p;
        double v_candidate, a_candidate, T;
        bool on_isochrone;
        vector<vector<double>> isochrone_right = {{v0, a0}};

        cout << "Start right branch of isochrone" << endl;
        while(true){
            p = next_point_right(isochrone_right);
            v_candidate = p[0];
            a_candidate = p[1];
            if(v_candidate < vmax + _step_width/2.){
                _isochrone.push_back({v_candidate, a_candidate});
                on_isochrone = false;
                while(!on_isochrone){
                    T = return_time_point_isochrone(v_candidate, a_candidate);
                    if (abs(T - _isi) / _isi > 0.001) {
                        if (T > _isi) {
                            cout << "ISI too long " << abs(T - _isi) / _isi << " " << "v: " << v_candidate << " a: "
                                 << a_candidate << endl;
                            a_candidate = a_candidate - max(0.001, 0.1 * abs(T - _isi) / _isi);
                        } else {
                            cout << "ISI too short " << abs(T - _isi) / _isi << " " << "v: " << v_candidate << " a: "
                                 << a_candidate << endl;
                            a_candidate = a_candidate + max(0.001, 0.1 * abs(T - _isi) / _isi);
                        }
                        if(a_candidate < 0.){

                        }
                        else{
                            _isochrone.back() = {v_candidate, a_candidate};
                        }
                    } else {
                        cout << "on point" << endl;
                        on_isochrone = true;
                        if (abs(v_candidate - vmax) < 0.001) {
                            _start_point_upper_branch = {0 , a_candidate + _deltaa};
                        }
                        isochrone_right.push_back({v_candidate, a_candidate});
                    }
                }
            }
            else{
                return isochrone_right;
            }
        }
    }
};


int main() {
    float mu = 2.0;
    float D = 0.5;
    float tau_a = 2.0;
    float delta_a = 1.0;
    double pi = 3.14159265;

    double dt = 0.0001;

    struct passwd *pw = getpwuid(getuid());
    const char *homedir = pw->pw_dir;

    double v0 = 0;
    double a0 = 2;

    Adaptive_leaky_integrate_and_fire alif(v0, a0, mu, D, tau_a, delta_a);

    double T = alif.period();
    vector<vector<double>> limit_cycle = alif.get_limit_cycle();
    double phase = 2.*pi/2.;
    int idx_phase = static_cast<int>(limit_cycle.size()*(phase/(2*pi)));
    vector<double> p = limit_cycle[idx_phase];
    double v_lc = p[0];
    double a_lc = p[1];
    cout << "v lc: " << v_lc << " a lc: " << a_lc << endl;
    double step_width = 0.05;

    Bob_the_isochrone_builder iso(alif, T, v_lc, a_lc, step_width);
    double vres = 0;
    double vmax = 1;
    double vmin = -5;

    // Construct central branch (that passes the limit cycle)
    vector<vector<double>> iso_left_b1 = iso.construct_isochrone_branch_left_from_point(v_lc, a_lc, vres, -0.01, 0);
    vector<vector<double>> iso_right_b1 = iso.construct_isochrone_branch_right_from_point(v_lc, a_lc, vmax);

    vector<vector<vector<double>>> branches;
    vector<vector<double>> branch;
    for(auto ele: iso_left_b1){
        branch.push_back(ele);
    }
    for(auto ele: iso_right_b1){
        branch.push_back(ele);
    }
    branches.push_back(branch);

    // Construct lower branches from vT to vR.
    for(int i = 0; i<2; i++){
        p = iso.get_start_point_lower_branch();
        v0 = p[0];
        a0 = p[1];
        if(a0 < 0) {
            break;
        }
        iso.insert_to_isochrone(0, v0, a0);
        vector<vector<double>> iso_lower_branch = iso.construct_isochrone_branch_left_from_point(v0, a0, vres, -0.01, 0);
        branches.insert(branches.begin() + 0, iso_lower_branch);
    }


    // Expand lower branches from vR to -5.
    for(int i = 0; i < branches.size(); i++) {
        vector<vector<double>> current_branch = branches[i];
        p = current_branch[0];
        v0 = p[0];
        a0 = p[1];
        if(a0 < 0){
            continue;
        }
        vector<vector<double>> isochrone = iso.get_isochrone();
        int idx = 0;
        for(auto & ele: isochrone){
            if(ele == p){
                break;
            }
            idx +=1;
        }
        vector<vector<double>> iso_lower_left_branch = iso.construct_isochrone_branch_left_from_point(v0, a0, vres, vmin, idx);

        // It may happen that for some point on the new branch the adaptation becomes strong in a sense that the voltage
        // change of the voltage becomes negative because a > v. This can lead to the trajectory circumventing the isochrone.
        // This can be prevent by the following check:

        // Take the last point on the isochrone that was just constructed
        p = iso_lower_left_branch[0];
        v0 = p[0];
        a0 = p[1];
        // Let it evole backwards for one period
        p = alif.integrate_backwards_for_T(v0, a0);
        // The result is a new point p that that lies on the upper branch. The upper branch, yet to be constructed, should
        // not include points left to p because they would circumvent the lower isochrone.
        if(p[0] > vmin){
            vmin = p[0];
        }
        cout << p[0] << endl;

    }


    // Construct upper branches from vR to -5 and from vR to vMax
    for(int i = 0; i<3; i++) {
        p = iso.get_start_point_upper_branch();
        v0 = p[0];
        a0 = p[1];
        p = iso.check_point_on_isochrone(v0, a0);
        vector<vector<double>> isochrone = iso.get_isochrone();
        int idx = 0;
        for(auto & ele: isochrone){
            if(ele == p){
                break;
            }
            idx +=1;
        }
        v0 = p[0];
        a0 = p[1];
        vector<vector<double>> iso_upper_left_branch = iso.construct_isochrone_branch_left_from_point(v0, a0, vres, vmin, idx);
        vector<vector<double>> iso_upper_right_branch = iso.construct_isochrone_branch_right_from_point(v0, a0, vmax);

        p = iso_upper_left_branch[0];
        v0 = p[0];
        a0 = p[1];
        p = alif.integrate_backwards_for_T(v0, a0);
        if(p[0] > vmin) {
            vmin = p[0];
        }
    }

    /*p = iso.get_start_point_upper_branch();
    v0 = p[0];
    a0 = p[1];
    p = iso.check_point_on_isochrone(v0, a0);
    vector<vector<double>> isochrone = iso.get_isochrone();
    int idx = 0;
    for(auto & ele: isochrone){
        if(ele == p){
            break;
        }
        idx +=1;
    }
    v0 = p[0];
    a0 = p[1];
    vector<vector<double>> iso_upper_left_branch = iso.construct_isochrone_branch_left_from_point(v0, a0, vres, -2, idx);
    vector<vector<double>> iso_upper_right_branch = iso.construct_isochrone_branch_right_from_point(v0, a0, vmax);
    */

    std::string path;
    path = "../out/";
    char parameters[200];
    std::sprintf(parameters, "mu%.1f_taua%.1f_delta%.1f_phase%.2f", mu, tau_a, delta_a, phase);
    std::string out_file_iso;
    out_file_iso = path + "isochrone_" + parameters + ".dat";
    std::ofstream file_iso;
    file_iso.open(out_file_iso);

    if (!file_iso.is_open()) {
        std::cout << "Could not open file at: " << out_file_iso<< std::endl;
        std::cout << "This is where I am: " << std::string(homedir) << std::endl;
        return 1;
    }
    vector<vector<double>> isochrone = iso.get_isochrone();
    for(auto& elem: isochrone){
        file_iso << elem[0] << " " <<  elem[1] << "\n";
    }
    return 0;
}
