#ifndef LOCAL
#define NDEBUG
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <cfloat>
#include <ctime>
#include <cassert>
#include <map>
#include <utility>
#include <set>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <sstream>
#include <complex>
#include <stack>
#include <queue>
#include <numeric>
#include <list>
#include <iomanip>
#include <fstream>
#include <bitset>

using namespace std;

#define foreach(it, c) for (__typeof__((c).begin()) it=(c).begin(); it != (c).end(); ++it)
template <typename T> void print_container(ostream& os, const T& c) { const char* _s = " "; if (!c.empty()) { __typeof__(c.begin()) last = --c.end(); foreach (it, c) { os << *it; if (it != last) os << _s; } } }
template <typename T> ostream& operator<<(ostream& os, const vector<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const set<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const multiset<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const deque<T>& c) { print_container(os, c); return os; }
template <typename T, typename U> ostream& operator<<(ostream& os, const map<T, U>& c) { print_container(os, c); return os; }
template <typename T, typename U> ostream& operator<<(ostream& os, const pair<T, U>& p) { os << "(" << p.first << ", " << p.second << ")"; return os; }

template <typename T> void print(T a, int n, const string& split = " ") { for (int i = 0; i < n; i++) { cerr << a[i]; if (i + 1 != n) cerr << split; } cerr << endl; }
template <typename T> void print2d(T a, int w, int h, int width = -1, int br = 0) { for (int i = 0; i < h; ++i) { for (int j = 0; j < w; ++j) { if (width != -1) cerr.width(width); cerr << a[i][j] << ' '; } cerr << endl; } while (br--) cerr << endl; }
template <typename T> void input(T& a, int n) { for (int i = 0; i < n; ++i) cin >> a[i]; }
#define dump(v) (cerr << #v << ": " << v << endl)
// #define dump(v)

#define rep(i, n) for (int i = 0; i < (int)(n); ++i)
#define erep(i, n) for (int i = 0; i <= (int)(n); ++i)
#define all(a) (a).begin(), (a).end()
#define rall(a) (a).rbegin(), (a).rend()
#define clr(a, x) memset(a, x, sizeof(a))
#define sz(a) ((int)(a).size())
#define mp(a, b) make_pair(a, b)
#define ten(n) ((long long)(1e##n))

template <typename T, typename U> void upmin(T& a, const U& b) { a = min<T>(a, b); }
template <typename T, typename U> void upmax(T& a, const U& b) { a = max<T>(a, b); }
template <typename T> void uniq(T& a) { sort(a.begin(), a.end()); a.erase(unique(a.begin(), a.end()), a.end()); }
template <class T> string to_s(const T& a) { ostringstream os; os << a; return os.str(); }
template <class T> T to_T(const string& s) { istringstream is(s); T res; is >> res; return res; }
bool in_rect(int x, int y, int w, int h) { return 0 <= x && x < w && 0 <= y && y < h; }

typedef long long ll;
typedef pair<int, int> pint;
typedef unsigned long long ull;

const int DX[] = { 0, 1, 0, -1 };
const int DY[] = { 1, 0, -1, 0 };

class Population
{
public:
    static int queryRegion(int x1, int y1, int x2, int y2)
    {
        cout << '?' << endl;
        cout << x1 << ' ' << y1 << ' ' << x2 << ' ' << y2 << endl;
        cout.flush();
        int res;
        cin >> res;
        return res;
    }
};


struct Rect
{
    Rect(int low_x, int low_y, int high_x, int high_y) :
        low_x(low_x), low_y(low_y), high_x(high_x), high_y(high_y)
    {
    }

    int low_x, low_y, high_x, high_y;
    int area = 0;
    int pop = -1;
};
class PopulationMapping
{
public:
    vector<string> mapPopulation(int max_percentage, vector <string> world_map, int total_population)
    {
        const int max_population = (ll)total_population * max_percentage / 100;
        const int h = world_map.size(), w = world_map[0].size();
        vector<vector<bool>> selected(h, vector<bool>(w));

        const int box_h = h / 5;
        const int box_w = w / 5;
        vector<Rect> rects;
        for (int ly = 0; ly < h; ly += box_h)
        {
            for (int lx = 0; lx < w; lx += box_w)
            {
                int hy = min(ly + box_h, h);
                int hx = min(lx + box_w, w);
                Rect rect(lx, ly, hx, hy);
                for (int y = ly; y < hy; ++y)
                    for (int x = lx; x < hx; ++x)
                        if (world_map[y][x] == 'X')
                            ++rect.area;
                rects.push_back(rect);
            }
        }
        sort(all(rects), [](const Rect& a, const Rect& b){ return a.area > b.area; });

        const int qs = 22;
        assert(rects.size() >= qs);
        rects.erase(rects.begin() + qs, rects.end());
        rep(i, qs)
            rects[i].pop = Population::queryRegion(rects[i].low_x, rects[i].low_y, rects[i].high_x - 1, rects[i].high_y - 1);

        int best_area = 0;
        vector<Rect> best;
        rep(mask, 1 << qs)
        {
            int area = 0;
            int pop = 0;
            rep(i, qs)
            {
                if (mask >> i & 1)
                {
                    pop += rects[i].pop;
                    area += rects[i].area;
                }
            }

            if (pop <= max_population && area > best_area)
            {
                best_area = area;
                best.clear();
                rep(i, qs)
                    if (mask >> i & 1)
                        best.push_back(rects[i]);
            }
        }

        vector<string> res(h, string(w, '.'));
        for (auto& rect : best)
        {
            for (int y = rect.low_y; y < rect.high_y; ++y)
                for (int x = rect.low_x; x < rect.high_x; ++x)
                    res[y][x] = 'X';
        }
        return res;
    }
};


#ifdef LOCAL
int main()
{
    int max_percentage;
    cin >> max_percentage;
    int height;
    cin >> height;
    vector<string> world_map(height);
    input(world_map, height);
    int total_population;
    cin >> total_population;

    vector<string> ret = PopulationMapping().mapPopulation(max_percentage, world_map, total_population);
    cout << ret.size() << endl;
    for (auto& s : ret)
        cout << s << endl;
    cout.flush();
}
#endif
