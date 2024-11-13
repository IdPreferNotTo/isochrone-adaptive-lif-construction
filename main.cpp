#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <unistd.h>
#include <pwd.h>
#include <cmath>

using namespace std;


bool intersect(double ax1, double ay1, double ax2, double ay2, double bx1, double by1, double bx2, double by2) {
    // segment_a = [[ax1, ay1], [ax2, ay2]], segment_b = [[bx1, by1], [bx2, by2]]
    double aa, ab, ba, bb, xIntersec;
    // Check for overlap: If the largest x coordinate of segment _a is smaller than the smallest x coordinate
    // of segment b then there can be no intersection. (same for y)
    if (max(ax1, ax2) < min(bx1, bx2)) {
        return false;
    }
    if (max(ay1, ay2) < min(by1, by2)) {
        return false;
    }

    // If the is _a mutual interval calculate the x coordinate of that intersection point and check if it is in the interval.
    // Calculate fa(_a) = aa*x + ba = y and fb(x) = ab*x + bb = y
    aa = (ay1 - ay2) / (ax1 - ax2);  // slope of segment _a
    ab = (by1 - by2) / (bx1 - bx2);  // slope of segment b
    ba = ay1 - aa * ax1;  // y intercept of segment _a
    bb = by1 - ab * bx1;  // y intercep of segment b
    xIntersec = (bb - ba) / (aa - ab);  // x coordinate of intersection point
    if (xIntersec < max(min(ax1, ax2), min(bx1, bx2)) or xIntersec > min(max(ax1, ax2), max(bx1, bx2))) {
        return false;
    } else {
        return true;
    }
}


class adaptiveLeakyIntegrateAndFire {
private:
    double _v, _a;
    float _mu, _taua, _deltaa;
    double _isi;

    double _dt = 0.0001;


public:
    adaptiveLeakyIntegrateAndFire(double v0, double a0, float i, float tauA, float deltaA){
        _v = v0;
        _a = a0;
        _mu = i;
        _taua = tauA;
        _deltaa = deltaA;
        _isi = period();

    }

    bool integrateDeterministic() {
        _v += (_mu - _v - _a) * _dt;
        _a += (-_a / _taua) * _dt;
        if (_v > 1.0) {
            _v = 0;
            _a += _deltaa;
            return true;
        } else {
            return false;
        }
    }

    vector<double> integrateDeterministic(double v, double a) const {
        v += (_mu - v - a) * _dt;
        a += (-a / _taua) * _dt;
        if (v > 1.0) {
            v -= 1;
            a += _deltaa;
        }
        return {v, a};

    }

    vector<vector<double>> getLimitCycle() {
        int spikeCount;
        bool fired;

        spikeCount = 0;
        // Allow the adaptive LIF model to evolve for some time so that we are close to the limit cycle
        while (spikeCount < 100) {
            fired = integrateDeterministic();
            if (fired) {
                spikeCount += 1;
            }
        }
        // Integrate one more time (from reset to threshold) and save as limit cycle
        std::vector<std::vector<double>> limitCycle;
        std::vector<double> p;
        while (true) {
            p = {_v, _a};
            limitCycle.push_back(p);
            fired = integrateDeterministic();
            if (fired) {
                return limitCycle;
            }
        }
    }

    double period() {
        double v, a, t;
        int countSpikes;
        bool fired;

        v = 0;
        a = 0;
        countSpikes = 0;
        while (countSpikes < 100) {
            fired = integrateDeterministic();
            if (fired) {
                countSpikes += 1;
            }
        }
        t = 0;
        countSpikes = 0;
        while (countSpikes < 100) {
            fired = integrateDeterministic();
            t += _dt;
            if (fired) {
                countSpikes += 1;
            }
        }
        return t / countSpikes;
    }

    vector<double> integrateBackwardsForT(double v0, double a0) const {
        double v = v0;
        double a = a0;
        double t = 0;
        while (t < _isi) {
            t += _dt;
            v -= (_mu - v - a) * _dt;
            a -= (-a / _taua) * _dt;
        }
        return {v, a};
    }

    double getDeltaa() const {
        return _deltaa;
    }

    double getDt() const {
        return _dt;
    }

    double getA() const {
        return _a;
    }

    double getV() const {
        return _v;
    }

    void setA(double a) {
        a = a;
    }

    void setV(double v) {
        v = v;
    }

};


class isochroneBuilder {
private:
    adaptiveLeakyIntegrateAndFire &model;
    double _isi;
    double _v0, _a0;
    double _dt, _deltaa;
    double _stepWidth;
    vector<double> startPointLowerBranch, startPointUpperBranch;
    vector<vector<double>> isochrone;

public:
    isochroneBuilder(adaptiveLeakyIntegrateAndFire &model, double isi, double v0, double a0, double stepWidth)
            : model(model) {
        _isi = isi;
        _v0 = v0;
        _a0 = a0;
        _dt = model.getDt();
        _deltaa = model.getDeltaa();
        _stepWidth = stepWidth;
        isochrone.push_back({v0, a0});
    }

    void insertToIsochrone(int idx, double v0, double a0) {
        isochrone.insert(isochrone.begin() + idx, {v0, a0});
    }

    vector<double> getStartPointLowerBranch() const {
        return startPointLowerBranch;
    }

    vector<double> getStartPointUpperBranch() const {
        return startPointUpperBranch;
    }

    vector<vector<double>> getIsochrone() {
        return isochrone;
    }

    vector<double> nextPointRight(vector<vector<double>> &list) const {
        vector<double> p1, p0;
        double v1, a1, v0, a0, dv, da, vCand, aCand, dif;
        int n = list.size();
        if (list.size() > 1) {
            p1 = list[n - 1];
            p0 = list[n - 2];
            v1 = p1[0];
            a1 = p1[1];
            v0 = p0[0];
            a0 = p0[1];
            dv = v1 - v0;
            da = a1 - a0;
            da = _stepWidth * (da / abs(dv));
            dv = _stepWidth * (dv / abs(dv));
            vCand = v1 + dv;
            aCand = a1 + da;
        } else {
            p1 = list[n - 1];
            v1 = p1[0];
            a1 = p1[1];
            dif = fmod(v1, 0.05);
            vCand = v1 + 0.05 - dif;
            aCand = a1;
        }
        return {vCand, aCand};
    }

    vector<double> nextPointLeft(vector<vector<double>> &list) const {
        vector<double> p1, p0;
        double v1, a1, v0, a0, dv, da, vCand, aCand, dif;
        if (list.size() > 1) {
            p1 = list[0];
            p0 = list[1];
            v1 = p1[0];
            a1 = p1[1];
            v0 = p0[0];
            a0 = p0[1];
            dv = v1 - v0;
            da = a1 - a0;
            da = _stepWidth * (da / abs(dv));
            dv = _stepWidth * (dv / abs(dv));
            vCand = v1 + dv;
            aCand = a1 + da;
        } else {
            p1 = list[0];
            v1 = p1[0];
            a1 = p1[1];
            dif = fmod(v1, 0.05);
            if (dif <= 0.0001) {
                dif = 0.05;
            }
            vCand = v1 - dif;
            aCand = a1;
        }
        return {vCand, aCand};
    }

    double returnTimePointIsochrone(double v, double a) {
        double t;
        vector<double> p, pBefore, pAfter;
        double vBefore, aBefore, vAfter, aAfter;
        double vIsoBefore, aIsoBefore, vIsoAfter, aIsoAfter;
        bool passedIsochrone, fired;
        int countRef;

        t = 0;
        countRef = 0;
        // Let the point (_v, _a) evolve until it hits the presumed isochrone
        while (true) {
            vBefore = v;
            aBefore = a;
            p = model.integrateDeterministic(v, a);
            v = p[0];
            a = p[1];
            if (abs(v - vBefore) > 0.5) {
                vAfter = v + 1.;
                aAfter = a - _deltaa;
            } else {
                vAfter = v;
                aAfter = a;
            }
            t += _dt;
            countRef += 1;
            // For every segment of the isochrone check if this segment and the segment that describes
            //  the change of _v, _a ((v_tmp, a_tmp), (_v, _a)) intersect.
            for (int i = 0; i < isochrone.size() - 1; i++) {
                pBefore = isochrone[i];
                pAfter = isochrone[i + 1];
                vIsoBefore = pBefore[0];
                aIsoBefore = pBefore[1];
                vIsoAfter = pAfter[0];
                aIsoAfter = pAfter[1];
                if (std::abs(vIsoAfter - vIsoBefore) > 0.5) {
                    continue;
                }
                passedIsochrone = intersect(vBefore, aBefore, vAfter, aAfter, vIsoBefore, aIsoBefore,
                                            vIsoAfter, aIsoAfter);
                if (passedIsochrone and countRef > 1) {
                    return t;
                }
            }
        }
    }

    vector<double> checkPointOnIsochrone(double vCandidate, double aCandidate) {
        cout << "Start new branch" << endl;
        bool onIsochrone;
        double t;

        isochrone.push_back({vCandidate, aCandidate});
        onIsochrone = false;
        while (!onIsochrone) {
            t = returnTimePointIsochrone(vCandidate, aCandidate);
            if (abs(t - _isi) / _isi > 0.002) {
                if (t > _isi) {
                    cout << std::setprecision(4) << "ISI too long " << abs(t - _isi) / _isi << " " << "_v: "
                         << vCandidate << " _a: "
                         << aCandidate << endl;
                    aCandidate = aCandidate - max(0.001, 0.1 * abs(t - _isi) / _isi);
                } else {
                    cout << std::setprecision(4) << "ISI too short " << abs(t - _isi) / _isi << " " << "_v: "
                         << vCandidate << " _a: "
                         << aCandidate << endl;
                    aCandidate = aCandidate + max(0.001, 0.1 * abs(t - _isi) / _isi);
                }
                isochrone.back() = {vCandidate, aCandidate};
            } else {
                cout << "on point" << endl;
                onIsochrone = true;
            }
        }
        return {vCandidate, aCandidate};
    }

    vector<vector<double>>
    constructIsochroneBranchLeftFromPoint(double v0, double a0, double vres, double vmin, int insertPoint) {
        // Add initial points (_v0, _a0) to left part of isochrone
        vector<double> p;
        double vCandidate, aCandidate, T;
        bool onIsochrone;
        vector<vector<double>> isochroneLeft = {{v0, a0}};

        cout << "Start left branch of isochrone" << endl;
        while (true) {
            p = nextPointLeft(isochroneLeft);
            vCandidate = p[0];
            aCandidate = p[1];
            double aBeforeOnIsochrone = isochroneLeft[0][1];

            if (vCandidate > vmin and aBeforeOnIsochrone > 0) {
                // Insert (_v0, _a0) to the general isochrone this is important because i calculate the return time to this
                // isochrone. Note that (_v0, _a0) is added at the beginning of isochrone.
                isochrone.insert(isochrone.begin() + insertPoint, {vCandidate, aCandidate});
                onIsochrone = false;
                // Calculate return time from (_v0, _a0) to isochrone and adjust (_v0, _a0) until |T - ISI| / ISI < 0.0005
                while (!onIsochrone) {
                    T = returnTimePointIsochrone(vCandidate, aCandidate);
                    if (abs(T - _isi) / _isi > 0.001) {
                        if (T > _isi) {
                            cout << std::setprecision(4) << "ISI too long " << abs(T - _isi) / _isi << " " << "_v: "
                                 << vCandidate << " _a: "
                                 << aCandidate << endl;
                            aCandidate = aCandidate - max(0.001, 0.1 * abs(T - _isi) / _isi);
                        } else {
                            cout << std::setprecision(4) << "ISI too short " << abs(T - _isi) / _isi << " " << "_v: "
                                 << vCandidate << " _a: "
                                 << aCandidate << endl;
                            aCandidate = aCandidate + max(0.001, 0.1 * abs(T - _isi) / _isi);
                        }
                        isochrone[insertPoint] = {vCandidate, aCandidate};
                    } else {
                        cout << "on point" << endl;
                        onIsochrone = true;
                        if (abs(vCandidate - vres) < 0.001) {
                            startPointLowerBranch = {1, aCandidate - _deltaa};
                        }
                        isochroneLeft.insert(isochroneLeft.begin() + 0, {vCandidate, aCandidate});
                    }
                }
            } else {
                return isochroneLeft;
            }
        }
    }

    vector<vector<double>> constructIsochroneBranchRightFromPoint(double v0, double a0, double vmax) {
        vector<double> p;
        double vCandidate, aCandidate, T;
        bool onIsochrone;
        vector<vector<double>> isochroneRight = {{v0, a0}};

        cout << "Start right branch of isochrone" << endl;
        while (true) {
            p = nextPointRight(isochroneRight);
            vCandidate = p[0];
            aCandidate = p[1];
            if (vCandidate < vmax + _stepWidth / 2.) {
                isochrone.push_back({vCandidate, aCandidate});
                onIsochrone = false;
                while (!onIsochrone) {
                    T = returnTimePointIsochrone(vCandidate, aCandidate);
                    if (abs(T - _isi) / _isi > 0.001) {
                        if (T > _isi) {
                            cout << "ISI too long " << abs(T - _isi) / _isi << " " << "_v: " << vCandidate << " _a: "
                                 << aCandidate << endl;
                            aCandidate = aCandidate - max(0.001, 0.1 * abs(T - _isi) / _isi);
                        } else {
                            cout << "ISI too short " << abs(T - _isi) / _isi << " " << "_v: " << vCandidate << " _a: "
                                 << aCandidate << endl;
                            aCandidate = aCandidate + max(0.001, 0.1 * abs(T - _isi) / _isi);
                        }
                        if (aCandidate < 0.) {

                        } else {
                            isochrone.back() = {vCandidate, aCandidate};
                        }
                    } else {
                        cout << "on point" << endl;
                        onIsochrone = true;
                        if (abs(vCandidate - vmax) < 0.001) {
                            startPointUpperBranch = {0, aCandidate + _deltaa};
                        }
                        isochroneRight.push_back({vCandidate, aCandidate});
                    }
                }
            } else {
                return isochroneRight;
            }
        }
    }
};


int main() {
    float mu = 2.0;
    float tauA = 2.0;
    float deltaA = 1.0;
    double pi = 3.14159265;

    double dt = 0.0001;

    struct passwd *pw = getpwuid(getuid());
    const char *homedir = pw->pw_dir;

    double v0 = 0;
    double a0 = 2;

    adaptiveLeakyIntegrateAndFire alif(v0, a0, mu, tauA, deltaA);

    double T = alif.period();
    vector<vector<double>> limitCycle = alif.getLimitCycle();
    double phase = 2. * pi / 2.;
    int idxPhase = static_cast<int>(limitCycle.size() * (phase / (2 * pi))) - 1.;
    cout << limitCycle.size() << endl;
    vector<double> p = limitCycle[idxPhase];
    double vLc = p[0];
    double aLc = p[1];
    cout << "_v lc: " << vLc << " _a lc: " << aLc << endl;
    double stepWidth = 0.05;

    isochroneBuilder iso(alif, T, vLc, aLc, stepWidth);
    double vres = 0;
    double vmax = 1;
    double vmin = -10;

    // Construct central branch (that passes the limit cycle)
    vector<vector<double>> isoLeftB1 = iso.constructIsochroneBranchLeftFromPoint(vLc, aLc, vres, -0.01, 0);
    vector<vector<double>> isoRightB1 = iso.constructIsochroneBranchRightFromPoint(vLc, aLc, vmax);

    vector<vector<vector<double>>> branches;
    vector<vector<double>> branch;
    for (auto ele: isoLeftB1) {
        branch.push_back(ele);
    }
    for (auto ele: isoRightB1) {
        branch.push_back(ele);
    }
    branches.push_back(branch);

    // Construct lower branches from vT to vR.
    for (int i = 0; i < 2; i++) {
        p = iso.getStartPointLowerBranch();
        v0 = p[0];
        a0 = p[1];
        if (a0 < 0) {
            break;
        }
        iso.insertToIsochrone(0, v0, a0);
        vector<vector<double>> isoLowerBranch = iso.constructIsochroneBranchLeftFromPoint(v0, a0, vres, -0.01,
                                                                                          0);
        branches.insert(branches.begin() + 0, isoLowerBranch);
    }


    // Expand lower branches from vR to -5.
    for (auto currentBranch: branches) {
        p = currentBranch[0];
        v0 = p[0];
        a0 = p[1];
        if (a0 < 0) {
            continue;
        }
        vector<vector<double>> isochrone = iso.getIsochrone();
        int idx = 0;
        for (auto &ele: isochrone) {
            if (ele == p) {
                break;
            }
            idx += 1;
        }
        vector<vector<double>> isoLowerLeftBranch = iso.constructIsochroneBranchLeftFromPoint(v0, a0, vres,
                                                                                              vmin, idx);

        // It may happen that for some point on the new branch the adaptation becomes strong in _a sense that the voltage
        // change of the voltage becomes negative because _a > _v. This can lead to the trajectory circumventing the isochrone.
        // This can be prevent by the following check:

        // Take the last point on the isochrone that was just constructed
        p = isoLowerLeftBranch[0];
        v0 = p[0];
        a0 = p[1];
        // Let it evole backwards for one period
        p = alif.integrateBackwardsForT(v0, a0);
        // The result is _a new point p that that lies on the upper branch. The upper branch, yet to be constructed, should
        // not include points left to p because they would circumvent the lower isochrone.
        if (p[0] > vmin) {
            vmin = p[0];
        }
        cout << p[0] << endl;

    }


    // Construct upper branches from vR to -5 and from vR to vMax
    for (int i = 0; i < 3; i++) {
        p = iso.getStartPointUpperBranch();
        v0 = p[0];
        a0 = p[1];
        p = iso.checkPointOnIsochrone(v0, a0);
        vector<vector<double>> isochrone = iso.getIsochrone();
        int idx = 0;
        for (auto &ele: isochrone) {
            if (ele == p) {
                break;
            }
            idx += 1;
        }
        v0 = p[0];
        a0 = p[1];
        vector<vector<double>> isoUpperLeftBranch = iso.constructIsochroneBranchLeftFromPoint(v0, a0, vres,
                                                                                              vmin, idx);
        vector<vector<double>> isoUpperRightBranch = iso.constructIsochroneBranchRightFromPoint(v0, a0, vmax);

        p = isoUpperLeftBranch[0];
        v0 = p[0];
        a0 = p[1];
        p = alif.integrateBackwardsForT(v0, a0);
        if (p[0] > vmin) {
            vmin = p[0];
        }
    }

    std::string path;
    path = "../out/";
    char parameters[200];
    std::sprintf(parameters, "_mu%.1f_taua%.1f_delta%.1f_phase%.2f", mu, tauA, deltaA, phase);
    std::string outFileIso;
    outFileIso = path + "isochrone_" + parameters + ".dat";
    std::ofstream fileIso;
    fileIso.open(outFileIso);

    if (!fileIso.is_open()) {
        std::cout << "Could not open file at: " << outFileIso << std::endl;
        std::cout << "This is where I am: " << std::string(homedir) << std::endl;
        return 1;
    }
    vector<vector<double>> isochrone = iso.getIsochrone();
    for (auto &elem: isochrone) {
        fileIso << elem[0] << " " << elem[1] << "\n";
    }
    return 0;
}
