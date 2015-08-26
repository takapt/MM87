#ifndef LOCAL
// #define NDEBUG
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


ull rdtsc()
{
#ifdef __amd64
    ull a, d;
    __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));
    return (d<<32) | a;
#else
    ull x;
    __asm__ volatile ("rdtsc" : "=A" (x));
    return x;
#endif
}
#ifdef LOCAL
const double CYCLES_PER_SEC = 3.30198e9;
#else
const double CYCLES_PER_SEC = 2.5e9;
#endif
double get_absolute_sec()
{
    return (double)rdtsc() / CYCLES_PER_SEC;
}
#ifdef _MSC_VER
#include <Windows.h>
    double get_ms() { return (double)GetTickCount64() / 1000; }
#else
#include <sys/time.h>
    double get_ms() { struct timeval t; gettimeofday(&t, NULL); return (double)t.tv_sec * 1000 + (double)t.tv_usec / 1000; }
#endif

#define USE_RDTSC
class Timer
{
private:
    double start_time;
    double elapsed;

#ifdef USE_RDTSC
    double get_sec() { return get_absolute_sec(); }
#else
    double get_sec() { return get_ms() / 1000; }
#endif

public:
    Timer() {}

    void start() { start_time = get_sec(); }
    double get_elapsed() { return elapsed = get_sec() - start_time; }
};

class Random
{
private:
    unsigned int  x, y, z, w;
public:
    Random(unsigned int x
             , unsigned int y
             , unsigned int z
             , unsigned int w)
        : x(x), y(y), z(z), w(w) { }
    Random()
        : x(123456789), y(362436069), z(521288629), w(88675123) { }
    Random(unsigned int seed)
        : x(123456789), y(362436069), z(521288629), w(seed) { }

    unsigned int next()
    {
        unsigned int t = x ^ (x << 11);
        x = y;
        y = z;
        z = w;
        return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    }

    int next_int() { return next(); }

    // [0, upper)
    int next_int(int upper) { return next() % upper; }

    // [low, high)
    int next_int(int low, int high) { return next_int(high - low) + low; }

    double next_double(double upper) { return upper * next() / UINT_MAX; }
    double next_double(double low, double high) { return next_double(high - low) + low; }

    template <typename T>
    int select(const vector<T>& ratio)
    {
        T sum = accumulate(ratio.begin(), ratio.end(), (T)0);
        double v = next_double(sum) + 1e-6;
        for (int i = 0; i < (int)ratio.size(); ++i)
        {
            v -= ratio[i];
            if (v <= 0)
                return i;
        }
        return (int)ratio.size() - 1;
    }
};
Random g_rand(3);


#ifdef LOCAL

#ifdef LOG_VIS_INPUT

ofstream log_input("input");
#define LOG_INPUT(s) {log_input << s << '\n'; log_input.flush();}

#else

#define LOG_INPUT(s)

#endif // LOG_VIS_INPUT
#endif // LOCAL

#ifdef LOCAL
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
        LOG_INPUT(res);
        return res;
    }
};
#endif

#ifdef LOCAL
int survey_size = 20;
#endif

int query_region(int low_x, int low_y, int high_x, int high_y)
{
    assert(low_x < high_x);
    assert(low_y < high_y);
    return Population::queryRegion(low_x, low_y, high_x - 1, high_y - 1);
}


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

class BitBoard
{
public:
    BitBoard(const vector<string>& s) :
        h(s.size()), w(s[0].size())
    {
        rep(y, h) rep(x, w)
            if (s[y][x] == 'X')
                f[index(x, y)] = true;
    }
    BitBoard(int w, int h) :
        h(h), w(w)
    {
    }

    bool get(int x, int y) const
    {
        return f[index(x, y)];
    }

    void set(int x, int y, bool val)
    {
        f[index(x, y)] = val;
    }

    int count(int low_x, int low_y, int high_x, int high_y) const
    {
        int c = 0;
        for (int y = low_y; y < high_y; ++y)
            for (int x = low_x; x < high_x; ++x)
                if (get(x, y))
                    ++c;
        return c;
    }

    int width() const { return w; }
    int height() const {return h; }

private:
    int index(int x, int y) const
    {
        assert(in_rect(x, y, w, h));
        return x | (y << 9);
    }

    const int h, w;
    bitset<512 * 512> f;
};


struct City
{
    City(int x, int y, int scale) :
        x(x), y(y), scale(scale)
    {
        assert(1 <= scale && scale <= 5);
    }

    int center_distance(int xx, int yy) const
    {
        return max(1, (abs(xx - x) + abs(yy - y)) / scale);
    }

    int x, y, scale;
};
namespace std
{
ostream& operator<<(ostream& os, const City& city)
{
    char buf[256];
    sprintf(buf, "(%3d, %3d, %d)", city.x, city.y, city.scale);
    os << buf;
    return os;
}
}

class PopMap
{
public:
    PopMap(int w, int h) :
        w(w), h(h), p{}
    {
    }
    PopMap(const BitBoard& world, const vector<City>& cities) :
        w(world.width()), h(world.height())
    {
        rep(y, h) rep(x, w)
        {
            if (world.get(x, y))
            {
                int center_dist = w + h;
                for (auto& city : cities)
                    upmin(center_dist, city.center_distance(x, y));
                int max_pop = (w + h) * 8 / center_dist;
                p[y][x] = max_pop / 2; // expect value
            }
            else
                p[y][x] = 0;
        }
    }
    PopMap(){}

    ll pop(int x, int y) const
    {
        assert(in_rect(x, y, w, h));
        return p[y][x];
    }

    ll pop(int low_x, int low_y, int high_x, int high_y) const
    {
        ll pp = 0;
        for (int y = low_y; y < high_y; ++y)
            for (int x = low_x; x < high_x; ++x)
                pp += pop(x, y);
        return pp;
    }

    void set_pop(int x, int y, int po)
    {
        assert(in_rect(x, y, w, h));
        p[y][x] = po;
    }

private:
    int w, h;
    int p[500][512];
};
#if 0
#define CORRECT_POP
#ifdef CORRECT_POP
PopMap correct_pop_map;
#endif
#endif

bool intersect(int a_low_x, int a_high_y, int b_low_x, int b_high_x)
{
    return a_low_x < b_high_x && b_low_x < a_high_y;
}

struct Region
{
    Region(int low_x, int low_y, int high_x, int high_y, int area, ll pop) :
        low_x(low_x), low_y(low_y), high_x(high_x), high_y(high_y), area(area), pop(pop)
    {
        assert(low_x < high_x);
        assert(low_y < high_y);
    }
    Region(int low_x, int low_y, int high_x, int high_y) :
        low_x(low_x), low_y(low_y), high_x(high_x), high_y(high_y), area(-1), pop(-1)
    {
        assert(low_x < high_x);
        assert(low_y < high_y);
    }
    Region(int low_x, int low_y, int high_x, int high_y, const BitBoard& world, const PopMap& pop_map) :
        low_x(low_x), low_y(low_y), high_x(high_x), high_y(high_y), area(world.count(low_x, low_y, high_x, high_y)), pop(pop_map.pop(low_x, low_y, high_x, high_y))
    {
        assert(low_x < high_x);
        assert(low_y < high_y);
    }
    Region() : Region(-114514, -1919, -114514 + 1, -1919 + 1){}

    int width() const { return high_x - low_x; }
    int height() const { return high_y - low_y; }

    Region moved(int dx, int dy, const BitBoard& world, const PopMap& pop_map) const
    {
        return Region(low_x + dx, low_y + dy, high_x + dx, high_y + dy, world, pop_map);
    }

    Region resized(int dx, int dy, const BitBoard& world, const PopMap& pop_map) const
    {
        int lx = low_x, hx = high_x, ly = low_y, hy = high_y;
        if (dx > 0)
            hx += dx;
        else if (dx < 0)
            lx += dx;

        if (dy > 0)
            hy += dy;
        else if (dy < 0)
            ly += dy;

        return Region(lx, ly, hx, hy, world, pop_map);
    }

    // true if either region contains the other
    bool intersect(const Region& other) const
    {
        return ::intersect(low_x, high_x, other.low_x, other.high_x) &&
               ::intersect(low_y, high_y, other.low_y, other.high_y);
    }

    // false if either region contains the other
    bool intersect_on_side(const Region& other) const
    {
        return intersect(other) &&
               !contain(other) && !other.contain(*this);
    }

    bool contain(const Region& other) const
    {
        return low_x <= other.low_x && other.high_x <= high_x &&
               low_y <= other.low_y && other.high_y <= high_y;
    }

    bool contain(int x, int y) const
    {
        return low_x <= x && x < high_x && low_y <= y && y < high_y;
    }

    bool operator==(const Region& other) const
    {
        return low_x == other.low_x &&
               low_y == other.low_y &&
               high_y == other.high_y &&
               high_y == other.high_y;
    }

    bool is_valid() const
    {
        if (low_x < 0 || low_y < 0 || low_x >= 500 || low_y >= 500)
            return false;
        if (area < 0 || pop < 0)
            return false;
        return true;
    }

    int low_x, low_y, high_x, high_y;
    int area;
    ll pop;
};
namespace std
{
ostream& operator<<(ostream& os, const Region& region)
{
    char buf[256];
    sprintf(buf, "[%3d, %3d), [%3d, %3d)", region.low_x, region.high_x, region.low_y, region.high_y);
    os << buf;
    return os;
}
}
int query_region(const Region& region)
{
    return query_region(region.low_x, region.low_y, region.high_x, region.high_y);
}

struct RegionNode
{
    RegionNode(const Region& region, RegionNode* parent, bool fixed = false) :
        region(region), parent(parent), fixed(fixed)
    {
    }

#ifndef NDEBUG
    ~RegionNode()
    {
        region = Region(-114514, -1919810, -114514 + 1, -1919810 + 1, -43062, -893931);
    }
#endif

    unique_ptr<RegionNode> copy_tree(RegionNode* parent = nullptr) const
    {
        unique_ptr<RegionNode> copy_node(new RegionNode(region, parent, fixed));
        copy_node->childs.reserve(childs.size());
        for (auto& child : childs)
            copy_node->childs.push_back(child->copy_tree(copy_node.get()));
        assert(copy_node);
        return copy_node;
    }

    const vector<unique_ptr<RegionNode>>& childs_nodes() const
    {
        return childs;
    }

    void add_child_region(const Region& reg, bool fixed = false)
    {
        assert(region.contain(reg));
        assert(!intersect_with_child(reg));
        childs.push_back(unique_ptr<RegionNode>(new RegionNode(reg, this, fixed)));
    }

    void remove_child_region(const RegionNode* node_to_remove)
    {
        auto it = childs.begin();
        while (it != childs.end() && it->get() != node_to_remove)
            ++it;
        assert(it != childs.end());
        childs.erase(it);
    }

    bool intersect_with_child(const Region& reg) const
    {
        for (auto& child : childs)
            if (reg.intersect(child->region))
                return true;
        return false;
    }

    bool intersect_with_child_on_side(const Region& reg) const
    {
        for (auto& child : childs)
            if (reg.intersect_on_side(child->region))
                return true;
        return false;
    }

    int area_except_childs() const
    {
        int area = region.area;
        for (auto& child : childs)
            area -= child->region.area;
        if (area < 0)
        {
            dump(region.area);
            dump(childs.size());
            dump(region);
            dump(childs[0]->region);
        }
        assert(area >= 0);
        return area;
    }

    ll pop_except_childs() const
    {
        int pop = region.pop;
        for (auto& child : childs)
            pop -= child->region.pop;
        if (pop < 0)
        {

            // FIXME
#if 0
#ifdef CORRECT_POP 
            dump(correct_pop_map.pop(region.low_x, region.low_y, region.high_x, region.high_y));
#endif
            cerr << "pop_except_childs" << endl;
            dump(pop);
#endif

#if 0
            dump(pop);
            dump(region.pop);
            dump(fixed);
            dump(parent->fixed);

            dump(parent->region);
            dump(region);

            dump(childs.size());
            for (auto& child : childs)
            {
                dump(child->region);
                dump(child->fixed);
            }
#endif
            upmax(pop, 0);
        }
        assert(pop >= 0);
        return pop;
    }

    void iterate(const function<void(const RegionNode*)>& callback) const
    {
        callback(this);
        for (auto& child : childs)
            ((const RegionNode*)child.get())->iterate(callback);
    }

    void iterate(const function<void(RegionNode*)>& callback)
    {
        callback(this);
        for (auto& child : childs)
            child->iterate(callback);
    }

    void iterate_pos_except_childs(const function<void(int, int)>& callback) const
    {
        for (int y = region.low_y; y < region.high_y; ++y)
        {
            for (int x = region.low_x; x < region.high_x; ++x)
            {
                bool contained_by_child = false;
                for (auto& child : childs)
                    contained_by_child |= child->region.contain(x, y);
                if (!contained_by_child)
                    callback(x, y);
            }
        }
    }

    const RegionNode* find_node_to_insert(const Region& region) const
    {
        const RegionNode* min_area_node = nullptr;
        int min_area = ten(9);
        iterate([&](const RegionNode* node) {
            if (node->region.area < min_area && node->region.contain(region))
            {
                min_area = node->region.area;
                min_area_node = node;
            }
        });
        return min_area_node;
    }

    RegionNode* find_node_to_insert(const Region& region)
    {
        RegionNode* min_area_node = nullptr;
        int min_area = ten(9);
        iterate([&](RegionNode* node) {
            if (node->region.area < min_area && node->region.contain(region))
            {
                min_area = node->region.area;
                min_area_node = node;
            }
        });
        return min_area_node;
    }

    int tree_size() const
    {
        int size = 0;
        iterate([&](const RegionNode*) { ++size; });
        return size;
    }

    vector<const RegionNode*> list_all(bool ignore_fixed) const
    {
        vector<const RegionNode*> nodes;
        iterate([&](const RegionNode* no) {
            if (!no->fixed || !ignore_fixed)
                nodes.push_back(no);
        });
        return nodes;
    }

    vector<RegionNode*> list_all(bool ignore_fixed)
    {
        vector<RegionNode*> nodes;
        iterate([&](RegionNode* no) {
            if (!no->fixed || !ignore_fixed)
                nodes.push_back(no);
        });
        return nodes;
    }

    vector<RegionNode*> list_leaves(bool ignore_fixed)
    {
        vector<RegionNode*> nodes;
        iterate([&](RegionNode* no) {
            if (no->childs.empty() && (!no->fixed || !ignore_fixed))
                nodes.push_back(no);
        });
        return nodes;
    }

    bool is_valid() const
    {
        if (intersect_with_child_on_side(region))
            return false;

        for (auto& child : childs)
        {
            if (!region.contain(child->region))
            {
                cerr << "Unko" << endl;
                return false;
            }
        }

        rep(j, childs.size()) rep(i, j)
            if (childs[i]->region.intersect(childs[j]->region))
                return false;

        return true;
    }

    bool is_valid_tree() const
    {
        bool valid = true;
        iterate([&](const RegionNode* no) {
            valid &= no->is_valid();
        });
        return valid;
    }

    Region region;
    RegionNode* parent;

    bool fixed = false;

private:
    vector<unique_ptr<RegionNode>> childs;
};

vector<City> search_cities(const BitBoard& world, const RegionNode* fixed_tree)
{
    assert(fixed_tree->tree_size() > 0);
    assert(fixed_tree->region.width() == world.width() && fixed_tree->region.height() == world.height());

    const ll total_population = fixed_tree->region.pop;

    const auto eval = [&](const PopMap& pop_map)
    {
        ll sum_diff = 0;
        const auto diff_pop = [&](const RegionNode* node)
        {
            ll pop = 0;
            node->iterate_pos_except_childs([&](int x, int y) {
                pop += pop_map.pop(x, y);
            });

            ll diff = abs(pop - node->pop_except_childs());
            sum_diff += diff;
        };
        fixed_tree->iterate(diff_pop);
        return sum_diff;
    };

    const int MAX_CITIES = 5 + max(0, (world.width() + world.height()) / 100 - 1);

    vector<City> best_cities;
    PopMap best_pop_map;
    ll best_score = ten(13);

    vector<City> current_cities;
    PopMap current_pop_map(world.width(), world.height());
    ll current_score = eval(current_pop_map);
    int last_best_i = -1;
    const int MAX_TRIES = 1000;
    rep(try_i, MAX_TRIES)
    {
        vector<City> next_cities = current_cities;
        const vector<int> next_ratio = {
            current_cities.size() > 0 ? 10 : 0, // move
            current_cities.size() > 0 ? 5 : 0, //change scale 
            current_cities.size() < MAX_CITIES ? 5 : 0, // generate
            current_cities.size() > 0 ? 5 : 0 // remove
        };
        switch (g_rand.select(next_ratio))
        {
            case 0:
            {
                // move
                assert(current_cities.size() > 0);
                int city_i = g_rand.next_int(next_cities.size());
                City& city = next_cities[city_i];

                const int MAX_MOVE = 20;
                int dx = g_rand.next_int(-min(MAX_MOVE, city.x), min(MAX_MOVE + 1, world.width() - city.x));
                int dy = g_rand.next_int(-min(MAX_MOVE, city.y), min(MAX_MOVE + 1, world.height() - city.y));
                city.x += dx;
                city.y += dy;
                assert(in_rect(city.x, city.y, world.width(), world.height()));
                if (!world.get(city.x, city.y))
                    continue;

                break;
            }
            case 1:
            {
                // change scale
                assert(current_cities.size() > 0);
                int city_i = g_rand.next_int(next_cities.size());
                City& city = next_cities[city_i];

                city.scale = g_rand.next_int(1, 5 + 1);

                break;
            }
            case 2:
            {
                // generate
                assert(current_cities.size() + 1 <= MAX_CITIES);

                int x, y;
                do {
                    x = g_rand.next_int(world.width());
                    y = g_rand.next_int(world.height());
                } while (!world.get(x, y));
                int scale = g_rand.next_int(1, 5 + 1);
                next_cities.push_back(City(x, y, scale));

                break;
            }
            case 3:
            {
                // remove
                assert(current_cities.size() > 0);
                int city_i = g_rand.next_int(next_cities.size());
                next_cities.erase(next_cities.begin() + city_i);

                break;
            }
            default:
                assert(false);
        }

        PopMap next_pop_map(world, next_cities);
        ll next_score = eval(next_pop_map);
        if (next_score < best_score)
        {
            best_score = next_score;
            best_cities = next_cities;
            best_pop_map = next_pop_map;

            last_best_i = try_i;

#ifdef CORRECT_POP
            fprintf(stderr, "%6d: %9lld, %5.4f\n", try_i, best_score, (double)best_score / total_population);
            ll truth_diff = 0;
            const int BOX_H = 20, BOX_W = 20;
            for (int ly = 0; ly < world.height(); ly += BOX_H)
            {
                for (int lx = 0; lx < world.width(); lx += BOX_W)
                {
                    ll predict_sum = 0;
                    ll truth_sum = 0;
                    for (int y = ly; y < min(ly + BOX_H, world.height()); ++y)
                    {
                        for (int x = lx; x < min(lx + BOX_W, world.width()); ++x)
                        {
                            if (world.get(x, y))
                            {
                                predict_sum += next_pop_map.pop(x, y);
                                truth_sum += correct_pop_map.pop(x, y);
                            }
                            else
                                assert(next_pop_map.pop(x, y) == 0);
                        }
                    }
                    ll diff = abs(predict_sum - truth_sum);
                    truth_diff += diff;
                }
            }
            fprintf(stderr, "%6d: %9lld, %5.4f\n", try_i, truth_diff, (double)truth_diff / total_population);
            cerr << endl;
//             dump(best_cities);
#endif
        }
        if (next_score < current_score)
        {
            current_score = next_score;
            current_cities = next_cities;
            current_pop_map = next_pop_map;
        }

        if (try_i - last_best_i > 0.2 * MAX_TRIES)
            break;
    }

#ifdef CORRECT_POP
    ll sum_diff = 0;
    rep(y, world.height()) rep(x, world.width())
    {
        if (world.get(x, y))
        {
            ll diff = abs(best_pop_map.pop(x, y) - correct_pop_map.pop(x, y));
            sum_diff += diff;
        }
    }
    const double diff_p = (double)sum_diff / total_population;
    assert(total_population == fixed_tree->region.pop);
    dump(total_population);
    dump(sum_diff);
    dump(diff_p);
    cerr << endl;
#endif

    return best_cities;
}


vector<ll> search_min_pop_for_area(const RegionNode* region_tree)
{
    assert(region_tree);
    const ll inf = ten(9);
    const int max_area = region_tree->region.area;
    vector<ll> min_pop(max_area + 1, inf); // min_pop[area] -> min pop
    min_pop[0] = 0;
    const auto callback = [&](const RegionNode* node)
    {
        const int area = node->area_except_childs();
        const ll pop = node->pop_except_childs();
        for (int i = max_area - area; i >= 0; --i)
            upmin(min_pop[i + area], min_pop[i] + pop);
    };
    region_tree->iterate(callback);
    return min_pop;
}
int search_max_area(const RegionNode* region_tree, const ll max_population)
{
    assert(region_tree);
    auto min_pop = search_min_pop_for_area(region_tree);
    for (int i = (int)min_pop.size() - 1; i >= 0; --i)
        if (min_pop[i] <= max_population)
            return i;
    return -1;
}

double search_max_score(const RegionNode* region_tree, const ll max_population)
{
    return search_max_area(region_tree, max_population) * pow(0.996, region_tree->tree_size() - 1);
}

BitBoard search_best_area_selection(const RegionNode* region_tree, const ll max_population)
{
    auto nodes = region_tree->list_all(false);
    const ll inf = ten(9);
    const int max_area = region_tree->region.area;
    vector<ll> min_pop(max_area + 1, inf); // min_pop[area] -> min pop
    vector<vector<int>> use_node(nodes.size() + 1, vector<int>(max_area + 1, -1));
    min_pop[0] = 0;
    rep(i, nodes.size())
    {
        const int area = nodes[i]->area_except_childs();
        const ll pop = nodes[i]->pop_except_childs();
        for (int a = max_area - area; a >= 0; --a)
        {
            if (min_pop[a] <= max_population)
            {
                use_node[i + 1][a] = use_node[i][a];

                ll npop = min_pop[a] + pop;
                if (npop < min_pop[a + area])
                {
                    min_pop[a + area] = npop;
                    use_node[i + 1][a + area] = i;
                }
            }
        }
    };
#ifndef NDEBUG
    auto mp = search_min_pop_for_area(region_tree);
    for (int i = 0; i <= max_area; ++i)
    {
        if (min_pop[i] <= max_population)
            assert(min_pop[i] == mp[i]);
    }
#endif

    int best_area = -1;
    for (int a = 0; a <= max_area; ++a)
        if (min_pop[a] <= max_population)
            best_area = a;
    assert(best_area != -1);
    dump(best_area);

    vector<bool> used_node(nodes.size());
    for (int i = use_node[nodes.size()][best_area], a = best_area; i >= 0; )
    {
        assert(0 <= i && i < nodes.size());
        assert(!used_node[i]);
        used_node[i] = true;
        a -= nodes[i]->area_except_childs();
        i = use_node[i][a];
    }

    // nodes is pre-order
    BitBoard selected(region_tree->region.high_x, region_tree->region.high_y);
    rep(i, nodes.size())
    {
        const bool u = used_node[i];
        const Region& region = nodes[i]->region;
        for (int y = region.low_y; y < region.high_y; ++y)
            for (int x = region.low_x; x < region.high_x; ++x)
                selected.set(x, y, u);
    }
    return selected;
}

unique_ptr<RegionNode> search_queries(const BitBoard& world, const PopMap& pop_map, const ll max_population, const RegionNode* fixed_region_tree)
{
#ifndef NDEBUG
    for (auto& node : fixed_region_tree->list_all(false))
    {
        assert(node->fixed);
        const Region& reg = node->region;
        assert(reg.area == world.count(reg.low_x, reg.low_y, reg.high_x, reg.high_y));
#ifdef CORRECT_POP
        assert(reg.pop == correct_pop_map.pop(reg.low_x, reg.low_y, reg.high_x, reg.high_y));
#endif
    }
#endif

    const int total_area = world.count(0, 0, world.width(), world.height());
    const ll total_population = pop_map.pop(0, 0, world.width(), world.height());
    int gens = 0;
    int gen_succ = 0;
    const auto generate_region = [&](const RegionNode* tree_) -> unique_ptr<RegionNode>
    {
        auto tree = tree_->copy_tree();
        auto nodes = tree->list_all(false);
        if (nodes.empty())
            return nullptr;
        ++gens;

        vector<ll> select_ratio(nodes.size());
        rep(i, nodes.size())
//             select_ratio[i] = nodes[i]->region.area;
            select_ratio[i] = (double)nodes[i]->region.area / total_area + (double)nodes[i]->pop_except_childs() / max_population;
        RegionNode* target_region_node = nodes[g_rand.select(select_ratio)];
        const Region& target = target_region_node->region;

        int low_x = g_rand.next_int(target.low_x, target.high_x);
        int low_y = g_rand.next_int(target.low_y, target.high_y);
        int high_x = g_rand.next_int(low_x, target.high_x) + 1;
        int high_y = g_rand.next_int(low_y, target.high_y) + 1;
        Region added_region(low_x, low_y, high_x, high_y);
        assert(target.contain(added_region));
        if (target_region_node->intersect_with_child(added_region))
            return nullptr;
        if (added_region == target)
            return nullptr;

        added_region.area = world.count(low_x, low_y, high_x, high_y);
        added_region.pop = pop_map.pop(low_x, low_y, high_x, high_y);
        target_region_node->add_child_region(added_region);

        assert(tree->is_valid_tree());

        ++gen_succ;
        return tree;
    };

    int rems = 0;
    int rem_succ = 0;
    const auto remove_region = [&](const RegionNode* tree_) -> unique_ptr<RegionNode>
    {
        ++rems;
        auto tree = tree_->copy_tree();
        auto nodes = tree->list_all(true);
        if (nodes.empty())
            return nullptr;

        RegionNode* remove_node = nodes[g_rand.next_int(nodes.size())];
        remove_node->parent->remove_child_region(remove_node);

        assert(tree->is_valid_tree());

        ++rem_succ;
        return tree;
    };

    int mvs = 0;
    int mv_succ = 0;
    const auto move_region = [&](const RegionNode* tree_) -> unique_ptr<RegionNode>
    {
        auto tree = tree_->copy_tree();
        auto nodes = tree->list_all(true);
        if (nodes.empty())
            return nullptr;
        ++mvs;

        RegionNode* move_region_node = nodes[g_rand.next_int(nodes.size())];
        const RegionNode* parent = move_region_node->parent;
        const Region& parent_region = parent->region;

        const Region& move_region = move_region_node->region;
        Region moved_region(-114514, -1919810, -114514 + 1, -1919810 + 1);
        const int MAX_MOVE = 20;
        switch (g_rand.next_int(4))
        {
            case 0: // right
                {
                    int max_dx = min(MAX_MOVE, parent_region.high_x - move_region.high_x);
                    for (auto& child : parent->childs_nodes())
                    {
                        if (child.get() != move_region_node)
                        {
                            const Region& reg = child->region;
                            if (reg.low_x >= move_region.high_x && intersect(move_region.low_y, move_region.high_y, reg.low_y, reg.high_y))
                                upmin(max_dx, reg.low_x - move_region.high_x);
                        }
                    }
                    for (auto& child : move_region_node->childs_nodes())
                        upmin(max_dx, child->region.low_x - move_region.low_x);

                    assert(max_dx >= 0);
                    if (max_dx == 0)
                        return nullptr;

                    int dx = 1 + g_rand.next_int(max_dx);
                    moved_region = move_region.moved(dx, 0, world, pop_map);
                    break;
                }
            case 1: // left
                {
                    int max_dx = min(MAX_MOVE, move_region.low_x - parent_region.low_x);
                    for (auto& child : parent->childs_nodes())
                    {
                        if (child.get() != move_region_node)
                        {
                            const Region& reg = child->region;
                            if (reg.high_x <= move_region.low_x && intersect(move_region.low_y, move_region.high_y, reg.low_y, reg.high_y))
                                upmin(max_dx, move_region.low_x - reg.high_x);
                        }
                    }
                    for (auto& child : move_region_node->childs_nodes())
                        upmin(max_dx, move_region.high_x - child->region.high_x);
                    assert(max_dx >= 0);
                    if (max_dx == 0)
                        return nullptr;

                    int dx = 1 + g_rand.next_int(max_dx);
                    moved_region = move_region.moved(-dx, 0, world, pop_map);
                    break;
                }
            case 2: // down
                {
                    int max_dy = min(MAX_MOVE, parent_region.high_y - move_region.high_y);
                    for (auto& child : parent->childs_nodes())
                    {
                        if (child.get() != move_region_node)
                        {
                            const Region& reg = child->region;
                            if (reg.low_y >= move_region.high_y && intersect(move_region.low_x, move_region.high_x, reg.low_x, reg.high_x))
                                upmin(max_dy, reg.low_y - move_region.high_y);
                        }
                    }
                    for (auto& child : move_region_node->childs_nodes())
                        upmin(max_dy, child->region.low_y - move_region.low_y);
                    assert(max_dy >= 0);
                    if (max_dy == 0)
                        return nullptr;

                    int dy = 1 + g_rand.next_int(max_dy);
                    moved_region = move_region.moved(0, dy, world, pop_map);
                    break;
                }
            case 3: // up
                {
                    int max_dy = min(MAX_MOVE, move_region.low_y - parent_region.low_y);
                    for (auto& child : parent->childs_nodes())
                    {
                        if (child.get() != move_region_node)
                        {
                            const Region& reg = child->region;
                            if (reg.high_y <= move_region.low_y && intersect(move_region.low_x, move_region.high_x, reg.low_x, reg.high_x))
                                upmin(max_dy, move_region.low_y - reg.high_y);
                        }
                    }
                    for (auto& child : move_region_node->childs_nodes())
                        upmin(max_dy, move_region.high_y - child->region.high_y);
                    assert(max_dy >= 0);
                    if (max_dy == 0)
                        return nullptr;

                    int dy = 1 + g_rand.next_int(max_dy);
                    moved_region = move_region.moved(0, -dy, world, pop_map);
                    break;
                }
            default:
                assert(false);
        }
        if (moved_region == parent_region)
            return nullptr;

        move_region_node->region = moved_region;

        assert(parent_region.contain(move_region));
        assert(tree->is_valid_tree());

        ++mv_succ;
        return tree;
    };

    int resizes = 0;
    int resize_succ = 0;
    const auto resize_region = [&](const RegionNode* tree_) -> unique_ptr<RegionNode>
    {
        auto tree = tree_->copy_tree();
        auto nodes = tree->list_all(true);
        if (nodes.empty())
            return nullptr;
        ++resizes;

        RegionNode* resize_region_node = nodes[g_rand.next_int(nodes.size())];
        const RegionNode* parent = resize_region_node->parent;
        const Region& parent_region = parent->region;

        const Region& resize_region = resize_region_node->region;
        Region resized_region(-114514, -1919810, -114514 + 1, -1919810 + 1);
        const int MAX_RESIZE = 20;
        switch (g_rand.next_int(4))
        {
            case 0: // right
                {
                    int max_dx = min(MAX_RESIZE, parent_region.high_x - resize_region.high_x);
                    for (auto& child : parent->childs_nodes())
                    {
                        if (child.get() != resize_region_node)
                        {
                            const Region& reg = child->region;
                            if (reg.low_x >= resize_region.high_x && intersect(resize_region.low_y, resize_region.high_y, reg.low_y, reg.high_y))
                                upmin(max_dx, reg.low_x - resize_region.high_x);
                        }
                    }
                    for (auto& child : resize_region_node->childs_nodes())
                        upmin(max_dx, child->region.low_x - resize_region.low_x);

                    assert(max_dx >= 0);
                    if (max_dx == 0)
                        return nullptr;

                    int dx = 1 + g_rand.next_int(max_dx);
                    resized_region = resize_region.resized(dx, 0, world, pop_map);
                    break;
                }
            case 1: // left
                {
                    int max_dx = min(MAX_RESIZE, resize_region.low_x - parent_region.low_x);
                    for (auto& child : parent->childs_nodes())
                    {
                        if (child.get() != resize_region_node)
                        {
                            const Region& reg = child->region;
                            if (reg.high_x <= resize_region.low_x && intersect(resize_region.low_y, resize_region.high_y, reg.low_y, reg.high_y))
                                upmin(max_dx, resize_region.low_x - reg.high_x);
                        }
                    }
                    for (auto& child : resize_region_node->childs_nodes())
                        upmin(max_dx, resize_region.high_x - child->region.high_x);
                    assert(max_dx >= 0);
                    if (max_dx == 0)
                        return nullptr;

                    int dx = 1 + g_rand.next_int(max_dx);
                    resized_region = resize_region.resized(-dx, 0, world, pop_map);
                    break;
                }
            case 2: // down
                {
                    int max_dy = min(MAX_RESIZE, parent_region.high_y - resize_region.high_y);
                    for (auto& child : parent->childs_nodes())
                    {
                        if (child.get() != resize_region_node)
                        {
                            const Region& reg = child->region;
                            if (reg.low_y >= resize_region.high_y && intersect(resize_region.low_x, resize_region.high_x, reg.low_x, reg.high_x))
                                upmin(max_dy, reg.low_y - resize_region.high_y);
                        }
                    }
                    for (auto& child : resize_region_node->childs_nodes())
                        upmin(max_dy, child->region.low_y - resize_region.low_y);
                    assert(max_dy >= 0);
                    if (max_dy == 0)
                        return nullptr;

                    int dy = 1 + g_rand.next_int(max_dy);
                    resized_region = resize_region.resized(0, dy, world, pop_map);
                    break;
                }
            case 3: // up
                {
                    int max_dy = min(MAX_RESIZE, resize_region.low_y - parent_region.low_y);
                    for (auto& child : parent->childs_nodes())
                    {
                        if (child.get() != resize_region_node)
                        {
                            const Region& reg = child->region;
                            if (reg.high_y <= resize_region.low_y && intersect(resize_region.low_x, resize_region.high_x, reg.low_x, reg.high_x))
                                upmin(max_dy, resize_region.low_y - reg.high_y);
                        }
                    }
                    for (auto& child : resize_region_node->childs_nodes())
                        upmin(max_dy, resize_region.high_y - child->region.high_y);
                    assert(max_dy >= 0);
                    if (max_dy == 0)
                        return nullptr;

                    int dy = 1 + g_rand.next_int(max_dy);
                    resized_region = resize_region.resized(0, -dy, world, pop_map);
                    break;
                }
            default:
                assert(false);
        }
        if (resized_region == parent_region)
            return nullptr;

        resize_region_node->region = resized_region;

        assert(parent_region.contain(resize_region));
        assert(tree->is_valid_tree());

        ++resize_succ;
        return tree;
    };


    const auto eval_official = [&](const RegionNode* tree)
    {
        auto min_pop = search_min_pop_for_area(tree);
        int max_area = 0;
        rep(i, min_pop.size())
            if (min_pop[i] <= max_population)
                max_area = i;

        return max_area * pow(0.996, tree->tree_size());
    };

    const auto calc_opt_area = [&]()
    {
        vector<int> pops;
        rep(y, world.height()) rep(x, world.width())
            if (world.get(x, y))
                pops.push_back(pop_map.pop(x, y));
        sort(all(pops));
        pops.push_back(max_population + 1);

        ll sum_pop = 0;
        rep(i, pops.size())
        {
            sum_pop += pops[i];
            if (sum_pop > max_population)
                return i;
        }
        assert(false);
    };
    const int opt_area = calc_opt_area();
    const auto eval = [&](const RegionNode* tree)
    {
        auto min_pop = search_min_pop_for_area(tree);
        int max_area = 0;
        double score = -1e9;
        rep(i, min_pop.size())
        {
//             if (min_pop[i] <= max_population * 1.2)
//                 max_area = i;
            double area_p = (double)i / opt_area;
            double pop_p = (double)min_pop[i] / max_population;
            double s = area_p - pop_p;
            upmax(score, s);
        }

        return score;// * pow(0.996, tree->tree_size());
    };

    unique_ptr<RegionNode> best_tree = fixed_region_tree->copy_tree();
    double best_score = eval_official(best_tree.get());
    rep(simu_i, 3)
    {
        int last_best_i = -1;
        unique_ptr<RegionNode> local_best_tree;
        double local_best_score = -1e9;

        unique_ptr<RegionNode> current_tree = fixed_region_tree->copy_tree();
        current_tree->iterate([](RegionNode* no) { no->fixed = true; });
        double current_score = eval(current_tree.get());

        const int MAX_TRIES = 6000;
        rep(try_i, MAX_TRIES)
        {
//             fprintf(stderr, "%6d: %3d, %f\n", try_i, current_tree->tree_size(), current_score);
            assert(current_tree);
            unique_ptr<RegionNode> next_tree;
            vector<int> next_ratio = { 10, 1, 50, 50 };
            if (current_tree->tree_size() >= 25)
                next_ratio[0] = 0;
            if (try_i > MAX_TRIES * 0.9)
                next_ratio = { 0, 1, 0, 0 };
            switch (g_rand.select(next_ratio))
            {
                case 0:
                    next_tree = move(generate_region(current_tree.get()));
                    break;
                case 1:
                    next_tree = move(remove_region(current_tree.get()));
                    break;
                case 2:
                    next_tree = move(resize_region(current_tree.get()));
                    break;
                case 3:
                    next_tree = move(move_region(current_tree.get()));
                    break;
                default:
                    assert(false);
            }
            if (next_tree)
            {
                double official_score = eval_official(next_tree.get());
                double score = eval(next_tree.get());
                if (official_score > local_best_score)
                {
                    local_best_score = official_score;
                    local_best_tree = next_tree->copy_tree();
                    last_best_i = try_i;
//                     fprintf(stderr, "%6d: %5.3f, %6d\n", try_i, local_best_score, search_max_area(local_best_tree.get(), max_population));
                }
                const double threshold = try_i < MAX_TRIES * 0.8 ? 0.9 : 1;
                if (score > current_score * threshold || official_score + 1e-9 > local_best_score)
                {
                    current_tree = move(next_tree);
                    current_score = score;
                }
            }
            const int TO_BEST = try_i < MAX_TRIES * 0.9 ? 200 : 5;
            if (local_best_tree && try_i - last_best_i > 200)
            {
                current_tree = local_best_tree->copy_tree();
                current_score = eval(local_best_tree.get());
                last_best_i = try_i;
            }
        }
        if (!local_best_tree)
            break;

        if (local_best_score > best_score)
        {
            best_score = local_best_score;
            best_tree = move(local_best_tree);
        }
    }

//     dump(gens);
//     dump(gen_succ);
//     dump(rems);
//     dump(rem_succ);
//     dump(mvs);
//     dump(mv_succ);
//     dump(resizes);
//     dump(resize_succ);
    return best_tree;
}

struct SearchRegionResult
{
    Region region;
    ll score = -1;

    bool is_valid() const
    {
        return score != -1;
    }
};
SearchRegionResult search_region_for_query(const RegionNode* tree, const RegionNode* prev_tree)
{
    SearchRegionResult result;
    assert(!result.is_valid());

    for (auto& node : tree->list_all(true))
    {
        const RegionNode* parent = prev_tree->find_node_to_insert(node->region);
        assert(parent != nullptr);

        ll score = min(node->region.area, parent->region.area - node->region.area);
        if (score > result.score)
        {
            result.score = score;
            result.region = node->region;
        }
    }
    return result;
}

BitBoard mark_predict_cities(const BitBoard& world, const vector<City>& cities)
{
    BitBoard selected(world.width(), world.height());
    for (auto& city : cities)
    {
        dump(city);
        const int s = city.scale * 2;
        int lx = max(0, city.x - s);
        int ly = max(0, city.y - s);
        int hx = min(world.width(), city.x + s);
        int hy = min(world.height(), city.y + s);
        query_region(lx, ly, hx, hy);
        for (int y = ly; y < hy; ++y)
            for (int x = lx; x < hx; ++x)
                selected.set(x, y, true);
    }
    return selected;
}

BitBoard select_area(const ll max_population, const BitBoard& world, const ll total_population)
{
    const int total_area = world.count(0, 0, world.width(), world.height());
    const Region whole_region(0, 0, world.width(), world.height(), total_area, total_population);
    unique_ptr<RegionNode> region_tree(new RegionNode(whole_region, nullptr, true));

//     while (true)
    vector<City> prepre;
    const int TRIES = 5;
    rep(_, TRIES)
    {
        Timer timer;

        timer.start();
        vector<City> predict_cities = search_cities(world, region_tree.get());
        fprintf(stderr, "search_cities time: %5.4f\n", timer.get_elapsed());
        prepre = predict_cities;

        PopMap predict_pop_map(world, predict_cities);

        assert(region_tree);
        timer.start();
        unique_ptr<RegionNode> searched_region_tree =
            search_queries(world, predict_pop_map, max_population, region_tree.get());
        fprintf(stderr, "search_queries time: %5.4f\n", timer.get_elapsed());
        assert(searched_region_tree);

        const double current_score = search_max_score(region_tree.get(), max_population);
        const double next_score = search_max_score(searched_region_tree.get(), max_population);
        if (next_score < current_score)
            break;

        SearchRegionResult result_for_query = search_region_for_query(searched_region_tree.get(), region_tree.get());
        if (!result_for_query.is_valid())
            break;

        Region region = result_for_query.region;
        assert(region.area == world.count(region.low_x, region.low_y, region.high_x, region.high_y));
        const ll pop = query_region(region);
        region.pop = pop;
#ifdef CORRECT_POP
        assert(region.pop == correct_pop_map.pop(region.low_x, region.low_y, region.high_x, region.high_y));
#endif

        RegionNode* node_to_insert = region_tree->find_node_to_insert(region);
        assert(node_to_insert != nullptr);
        node_to_insert->add_child_region(region, true);
    }
//     return mark_predict_cities(world, prepre);

//     for (auto& node : searched_region_tree->list_all(false))
//     {
//         const Region& region = node->region;
//         if (region.width() == world.width() && region.height() == world.height())
//             continue;
//
//         query_region(region.low_x, region.low_y, region.high_x, region.high_y);
//         dump(region);
// //         for (int y = region.low_y; y < region.high_y; ++y)
// //             for (int x = region.low_x; x < region.high_x; ++x)
// //                 selected.set(x, y, true);
//     }

    BitBoard selected = search_best_area_selection(region_tree.get(), max_population);
#ifdef CORRECT_POP
    int selected_pop = 0;
    rep(y, world.height()) rep(x, world.width())
        if (selected.get(x, y) && world.get(x, y))
            selected_pop += correct_pop_map.pop(x, y);
    dump(selected_pop);
    fprintf(stderr, "pop: %.2f\n", (double)selected_pop / total_population);
#endif
    return selected;
}

class PopulationMapping
{
public:
    vector<string> mapPopulation(int max_percentage, vector <string> world_map, int total_population)
    {
//         return simple_search(max_percentage, world_map, total_population);

        const ll max_population = (ll)max_percentage * total_population / 100;
        BitBoard world(world_map);
        const ll total_area = world.count(0, 0, world.width(), world.height());


        BitBoard selected = select_area(max_population, world, total_population);
        return make_string(selected);
//         BitBoard selected(world.width(), world.height());

        unique_ptr<RegionNode> root(new RegionNode(Region(0, 0, world.width(), world.height(), total_area, total_population), nullptr, true));
        auto predict_cities = search_cities(world, root.get());

        for (auto& city : predict_cities)
        {
            int lx = max(0, city.x - 3);
            int ly = max(0, city.y - 3);
            int hx = min(world.width(), city.x + 3);
            int hy = min(world.height(), city.y + 3);
            query_region(lx, ly, hx, hy);
            for (int y = ly; y < hy; ++y)
                for (int x = lx; x < hx; ++x)
                    selected.set(x, y, true);
        }

        return make_string(selected);
    }

private:
    vector<string> make_string(const BitBoard& selected)
    {
        vector<string> res(selected.height(), string(selected.width(), '.'));
        rep(y, selected.height()) rep(x, selected.width())
            if (selected.get(x, y))
                res[y][x] = 'X';
        return res;
    }
#ifdef LOCAL
    struct SurveyResult
    {
        int opt_area;
        vector<Rect> used_rects;
        vector<string> selected;
    };
    SurveyResult survey(int max_percentage, vector <string> world_map, int total_population)
    {
        const int max_population = (ll)total_population * max_percentage / 100;
        const int h = world_map.size(), w = world_map[0].size();
        BitBoard world(world_map);

        vector<Rect> smalls;
        const int S = survey_size;
        for (int y = 0; y < h; y += S)
        {
            for (int x = 0; x < w; x += S)
            {
                Rect r(x, y, min(x + S, w), min(y + S, h));
                r.area = world.count(x, y, r.high_x, r.high_y);
                if (r.area > 0)
                {
                    int p = query_region(x, y, r.high_x, r.high_y);
                    if (p)
                    {
                        r.pop = p;
                        assert(r.area > 0);
                        smalls.push_back(r);
                    }
                }
            }
        }
        sort(all(smalls), [](const Rect& a, const Rect& b){ return a.pop > b.pop; });
        ll ss = accumulate(all(smalls), 0, [](int s, const Rect& r){ return s + r.pop; });
        ll sum = 0;
        //         rep(i, smalls.size())
        //         {
        //             if (i % (smalls.size() / 20) == 0)
        //             {
        //                 fprintf(stderr, "%5d (%4.2f): %9lld (%4.2f)\n", i, (double)i / smalls.size(), sum, (double)sum / ss);
        //             }
        //             sum += smalls[i].pop;
        //         }

        const int total_area = world.count(0, 0, world.width(), world.height());

        Timer timer;
        timer.start();
        dump(smalls.size());
        int opt_area;
        int cur = 0, next = 1;
        int dp[2][512 * 512];
        vector<vector<short>> use(smalls.size() + 1, vector<short>(total_area + 1, -1));
        erep(j, total_area)
            dp[0][j] = max_population + 1;
        dp[0][0] = 0;
        rep(i, smalls.size())
        {
            erep(j, total_area)
                dp[next][j] = max_population + 1;

            const auto& r = smalls[i];
            for (int j = total_area; j >= 0; --j)
            {
                if (dp[cur][j] < dp[next][j])
                {
                    dp[next][j] = dp[cur][j];
                    use[i + 1][j] = use[i][j];
                }

                if (dp[cur][j] <= max_population && j + r.area <= total_area && dp[next][j + r.area] > dp[cur][j] + r.pop)
                {
                    upmin(dp[next][j + r.area], dp[cur][j] + r.pop);
                    use[i + 1][j + r.area] = i;
                }
            }

            swap(cur, next);
        }
        erep(j, total_area)
            if (dp[cur][j] <= max_population)
                opt_area = j;
        dump(timer.get_elapsed());
        vector<Rect> use_r;
        vector<bool> used_i(smalls.size());
        for (int i = use[smalls.size()][opt_area], j = opt_area; j > 0; )
        {
            assert(0 <= i && i < smalls.size());
            assert(!used_i[i]);
            used_i[i] = true;
            use_r.push_back(smalls[i]);
            assert(smalls[i].area > 0);
            j -= smalls[i].area;
            i = use[i][j];
        }
        assert(use_r.size() > 0);
        vector<string> res(h, string(w, '.'));
        for (auto& rect : use_r)
        {
            for (int y = rect.low_y; y < rect.high_y; ++y)
                for (int x = rect.low_x; x < rect.high_x; ++x)
                {
                    assert(res[y][x] == '.');
                    res[y][x] = 'X';
                }
        }
        return SurveyResult {
            .opt_area = opt_area,
            .used_rects = use_r,
            .selected = res
        };
    }


    vector<string> simple_search(int max_percentage, vector <string> world_map, int total_population)
    {
        const int max_population = (ll)total_population * max_percentage / 100;
        const int h = world_map.size(), w = world_map[0].size();
        vector<vector<bool>> selected(h, vector<bool>(w));

        BitBoard world(world_map);

        const SurveyResult survey_result = survey(max_percentage, world_map, total_population);
        int opt_area = survey_result.opt_area;
//         return survey_result.selected;

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
                rect.area = world.count(lx, ly, hx, hy);
                rects.push_back(rect);
            }
        }
        sort(all(rects), [](const Rect& a, const Rect& b){ return a.area > b.area; });

        const int qs = 22;
        assert(rects.size() >= qs);
        rects.erase(rects.begin() + qs, rects.end());
        rep(i, qs)
            rects[i].pop = query_region(rects[i].low_x, rects[i].low_y, rects[i].high_x, rects[i].high_y);

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

        fprintf(stderr, "my / opt: %8d %8d, %2.3f\n", best_area, opt_area, (double)best_area / opt_area);

        vector<string> res(h, string(w, '.'));
        for (auto& rect : best)
        {
            for (int y = rect.low_y; y < rect.high_y; ++y)
                for (int x = rect.low_x; x < rect.high_x; ++x)
                    res[y][x] = 'X';
        }
        return res;
    }
#endif
};


#ifdef LOCAL
int main(int argc, char* argv[])
{
    if (argc > 1)
        survey_size = to_T<int>(argv[1]);

    int max_percentage;
    cin >> max_percentage;
    LOG_INPUT(max_percentage);
    int height;
    cin >> height;
    LOG_INPUT(height);
    vector<string> world_map(height);
    input(world_map, height);
    for (auto& s : world_map)
        LOG_INPUT(s);
    int total_population;
    cin >> total_population;
    LOG_INPUT(total_population);

#ifdef CORRECT_POP
    int h, w;
    cin >> h >> w;
    LOG_INPUT(h);
    LOG_INPUT(w);
    correct_pop_map = PopMap(w, h);
    rep(y, h) rep(x, w)
    {
        int p;
        cin >> p;
        LOG_INPUT(p);
        correct_pop_map.set_pop(x, y, p);
    }
    assert(correct_pop_map.pop(0, 0, w, h) == total_population);
#endif

    vector<string> ret = PopulationMapping().mapPopulation(max_percentage, world_map, total_population);
    cout << ret.size() << endl;
    for (auto& s : ret)
        cout << s << endl;
    cout.flush();
}
#endif
