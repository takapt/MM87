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
Timer g_timer;
#ifdef LOCAL
const double G_TL_SEC = 40;
#else
const double G_TL_SEC = 20;
#endif

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

    // dirty hack
    template <typename T>
    int count(const T& region) const
    {
        return count(region.low_x, region.low_y, region.high_x, region.high_y);
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

    // dirty hack
    template <typename T>
    int pop(const T& region) const
    {
        return pop(region.low_x, region.low_y, region.high_x, region.high_y);
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
#define PREDICT_PERFECT
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
    int rect_area() const { return width() * height(); }

    pair<Region, Region> split_vertical(int split_x) const
    {
        return make_pair(Region(low_x, low_y, split_x, high_y), Region(split_x, low_y, high_x, high_y));
    }

    pair<Region, Region> split_horizontal(int split_y) const
    {
        return make_pair(Region(low_x, low_y, high_x, split_y), Region(low_x, split_y, high_x, high_y));
    }

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
               high_x == other.high_x &&
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

    unique_ptr<RegionNode> copy_fixed_tree(RegionNode* parent = nullptr) const
    {
        if (!fixed)
            return nullptr;

        unique_ptr<RegionNode> copy_node(new RegionNode(region, parent, fixed));
        for (auto& child : childs)
        {
            auto c = child->copy_fixed_tree(copy_node.get());
            if (c)
                copy_node->childs.push_back(move(c));
        }
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
            if (node->region.rect_area() < min_area && node->region.contain(region))
            {
                min_area = node->region.rect_area();
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
            if (node->region.rect_area() < min_area && node->region.contain(region))
            {
                min_area = node->region.rect_area();
                min_area_node = node;
            }
        });
        return min_area_node;
    }

    RegionNode* find_equal_region_node(const Region& region)
    {
        if (region == this->region)
            return this;
        for (auto& child : childs)
        {
            auto res = child->find_equal_region_node(region);
            if (res)
                return res;
        }
        return nullptr;
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

    vector<const RegionNode*> list_leaves(bool ignore_fixed) const
    {
        vector<const RegionNode*> nodes;
        iterate([&](const RegionNode* no) {
            if (no->childs.empty() && (!no->fixed || !ignore_fixed))
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

    vector<unique_ptr<RegionNode>> childs;
};

vector<City> search_cities(const BitBoard& world, const RegionNode* fixed_tree)
{
    assert(fixed_tree->tree_size() > 0);
    assert(fixed_tree->region.width() == world.width() && fixed_tree->region.height() == world.height());

    const ll total_population = fixed_tree->region.pop;

    const auto eval = [&](const PopMap& pop_map)//, vector<pair<Region, ll>>& diff_info)
    {
        ll sum_diff = 0;
        for (auto& node : fixed_tree->list_leaves(false))
        {
            if (fixed_tree->fixed)
            {
                ll pop = pop_map.pop(node->region);

                ll diff = abs(pop - node->region.pop);
                sum_diff += diff;

//                 diff_info.push_back(make_pair(node->region, diff));
            }
        };
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
    const int MAX_TRIES = 200;
    rep(try_i, MAX_TRIES)
    {
        vector<City> next_cities = current_cities;
        const vector<int> next_ratio = {
            current_cities.size() >= 5 ? 10 : 0, // move
            current_cities.size() >= 5 ? 5 : 0, //change scale 
            current_cities.size() < MAX_CITIES ? 5 : 0, // generate
            current_cities.size() > 5 ? 2 : 0 // remove
        };
        switch (g_rand.select(next_ratio))
        {
            case 0:
            {
                // move
                assert(current_cities.size() >= 5);
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
                assert(current_cities.size() >= 5);
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
                assert(current_cities.size() > 5);
                int city_i = g_rand.next_int(next_cities.size());
                next_cities.erase(next_cities.begin() + city_i);

                break;
            }
            default:
                assert(false);
        }
        if (current_cities.size() < 5)
        {
            current_score = 1e9;
            current_cities = next_cities;
            continue;
        }

        PopMap next_pop_map(world, next_cities);
        ll next_score = eval(next_pop_map);
        if (next_cities.size() >= 5 && next_score < best_score)
        {
            best_score = next_score;
            best_cities = next_cities;
            best_pop_map = next_pop_map;

            last_best_i = try_i;

#if 0
            fprintf(stderr, "%6d: %9lld, %5.4f\n", try_i, best_score, (double)best_score / total_population);
#ifdef CORRECT_POP
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
#endif
        }
        bool prob_next = (next_score - current_score) * ((double)try_i / MAX_TRIES) < total_population * 0.2;
//         prob_next = false;
        if (next_score < current_score || prob_next)
        {
            current_score = next_score;
            current_cities = next_cities;
            current_pop_map = next_pop_map;
        }

        if (best_cities.size() >= 5 && try_i - last_best_i > 0.2 * MAX_TRIES)
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
//     dump(total_population);
//     dump(sum_diff);
//     dump(diff_p);
//     cerr << endl;
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

struct SearchBestRegionSelectionResult
{
    SearchBestRegionSelectionResult(const vector<const RegionNode*>& nodes, const vector<bool>& selected_node, ll max_area) :
        nodes(nodes), selected_node(selected_node), max_area(max_area)
    {
        assert(nodes.size() == selected_node.size());
    }

    vector<const RegionNode*> selected_nodes() const
    {
        vector<const RegionNode*> res;
        rep(i, nodes.size())
            if (selected_node[i])
                res.push_back(nodes[i]);
        return res;
    }

    vector<const RegionNode*> not_selected_nodes() const
    {
        vector<const RegionNode*> res;
        rep(i, nodes.size())
            if (selected_node[i])
                res.push_back(nodes[i]);
        return res;
    }

    vector<const RegionNode*> nodes;
    vector<bool> selected_node;
    ll max_area;
};
SearchBestRegionSelectionResult search_best_region_selection(const RegionNode* region_tree, const ll max_population)
{
    auto nodes = region_tree->list_all(false);
    const ll inf = ten(9);
    const int max_area = region_tree->region.area;
    vector<ll> min_pop(max_area + 1, inf); // min_pop[area] -> min pop
    vector<vector<int>> prev_selected_node(nodes.size() + 1, vector<int>(max_area + 1, -1));
    min_pop[0] = 0;
    rep(i, nodes.size())
    {
        const int area = nodes[i]->area_except_childs();
        const ll pop = nodes[i]->pop_except_childs();
        for (int a = max_area - area; a >= 0; --a)
        {
            if (min_pop[a] <= max_population)
            {
                prev_selected_node[i + 1][a] = prev_selected_node[i][a];

                ll npop = min_pop[a] + pop;
                if (npop < min_pop[a + area])
                {
                    min_pop[a + area] = npop;
                    prev_selected_node[i + 1][a + area] = i;
                }
            }
        }
    };

    int best_area = -1;
    for (int a = 0; a <= max_area; ++a)
        if (min_pop[a] <= max_population)
            best_area = a;
    assert(best_area != -1);
//     dump(best_area);

    vector<bool> selected_node(nodes.size());
    for (int i = prev_selected_node[nodes.size()][best_area], a = best_area; i >= 0; )
    {
        assert(0 <= i && i < nodes.size());
        assert(!selected_node[i]);
        selected_node[i] = true;
        a -= nodes[i]->area_except_childs();
        i = prev_selected_node[i][a];
    }

    return SearchBestRegionSelectionResult(move(nodes), move(selected_node), best_area);
}

BitBoard search_best_area_selection(const RegionNode* region_tree, const ll max_population)
{
    auto result = search_best_region_selection(region_tree, max_population);

    // nodes is pre-order
    BitBoard selected(region_tree->region.high_x, region_tree->region.high_y);
    rep(i, result.nodes.size())
    {
        const bool u = result.selected_node[i];
        const Region& region = result.nodes[i]->region;
        for (int y = region.low_y; y < region.high_y; ++y)
            for (int x = region.low_x; x < region.high_x; ++x)
                selected.set(x, y, u);
    }
    return selected;
}

struct SearchRegionResult
{
    Region region;
    double score = -114514;

    bool is_valid() const
    {
        return score > -1;
    }
};

SearchRegionResult search_region_for_query(const BitBoard& world, const ll max_population, const RegionNode* current_tree, const function<double(const Region& region)>& expect_pop)
{
    const ll total_area = current_tree->region.area;
    const ll total_population = current_tree->region.pop;

//     const auto expect_pop = [&](const Region& region)
//     {
//         double sum_pop = 0;
//         for (auto& pop_map : predict_pop_maps)
//             sum_pop += pop_map.pop(region);
//         return sum_pop / predict_pop_maps.size();
//     };

    auto tree = current_tree->copy_tree();
    for (auto& node : tree->list_all(true))
        node->region.pop = expect_pop(node->region);

    SearchBestRegionSelectionResult region_selection_result = search_best_region_selection(tree.get(), max_population);
    auto selected_nodes = region_selection_result.selected_nodes();
    auto not_selected_nodes = region_selection_result.not_selected_nodes();
    const auto count_need_queries = [](const vector<const RegionNode*>& nodes)
    {
        int c = nodes.size();
        for (auto& node : nodes)
            if (node->fixed)
                --c;
        return c;
    };
    const int selected_need_queries = count_need_queries(selected_nodes);
    const int not_selected_need_queries = count_need_queries(not_selected_nodes);

    SearchRegionResult result;
    if (selected_need_queries == 0 || not_selected_need_queries == 0)
        return result;

    auto& cand_nodes = selected_need_queries <= not_selected_need_queries ? selected_nodes : not_selected_nodes;
    const auto eval = [&](const RegionNode* node)
    {
        double area_score = (double)node->region.area / total_area;
        double pop_score = (double)node->region.pop / total_population;
        double score = area_score + pop_score;
        return score;
    };
    for (auto& node : cand_nodes)
    {
        if (node->fixed)
            continue;

        double score = eval(node);
        if (score > result.score)
        {
            result.score = score;
            result.region = node->region;
        }
    }
    return result;
}

unique_ptr<RegionNode> build_div_reiong_tree(const BitBoard& world, const ll total_population)
{
    const int total_area = world.count(0, 0, world.width(), world.height());
    const Region whole_region(0, 0, world.width(), world.height(), total_area, total_population);
    unique_ptr<RegionNode> region_tree(new RegionNode(whole_region, nullptr, true));

    int DIV_SIZE = 20;
    while (DIV_SIZE * DIV_SIZE * 15 < total_area)
        ++DIV_SIZE;
    dump(DIV_SIZE);
    const int DIV_H = DIV_SIZE;
    const int DIV_W = DIV_SIZE;
    for (int ly = 0; ly < world.height(); ly += DIV_H)
    {
        for (int lx = 0; lx < world.width(); lx += DIV_W)
        {
            int hy = min(world.height(), ly + DIV_H);
            int hx = min(world.width(), lx + DIV_W);
            Region region(lx, ly, hx, hy);
            region.area = world.count(region);
            region_tree->add_child_region(region);
        }
    }
    return region_tree;
}

BitBoard mark_predict_cities(const BitBoard& world, const vector<City>& cities)
{
    BitBoard selected(world.width(), world.height());
    for (auto& city : cities)
    {
        dump(city);
        const int s = city.scale * 3;
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
    unique_ptr<RegionNode> region_tree = build_div_reiong_tree(world, total_population);

    vector<City> prepre;
    const int MAX_QUERIES = 40;
    rep(query_i, MAX_QUERIES)
    {
        if (g_timer.get_elapsed() > G_TL_SEC * 0.9)
            break;

        vector<vector<City>> predict_cities;
        vector<PopMap> predict_pop_maps;
#ifdef PREDICT_PERFECT
        const int MAX_PREDICTS = 1;
        predict_pop_maps = {correct_pop_map};
#else
        const int MAX_PREDICTS = 1;
        rep(i, MAX_PREDICTS)
        {
            auto cities = search_cities(world, region_tree.get());
            predict_cities.push_back(cities);
            predict_pop_maps.push_back(PopMap(world, cities));
        }
        if (query_i == 35)
        {
            assert(predict_cities.size() >= 1);
            return mark_predict_cities(world, predict_cities[0]);
        }
#endif

        const auto expect_pop = [&](const Region& region)
        {
            double sum_pop = 0;
            for (auto& pop_map : predict_pop_maps)
                sum_pop += pop_map.pop(region);
            return sum_pop / predict_pop_maps.size();
        };

        SearchRegionResult result_for_query = search_region_for_query(world, max_population, region_tree.get(), expect_pop);
        if (!result_for_query.is_valid())
            break;

        const Region& region = result_for_query.region;

        auto fixed_tree = region_tree->copy_fixed_tree();
        const double current_score = search_max_area(fixed_tree.get(), max_population) * pow(0.996, query_i);
        // TODO: max expect_score if use n queries. now 1 query.
        auto next_tree = region_tree->copy_tree();
        auto node = next_tree->find_equal_region_node(region);
        node->fixed = true;
        node->region.pop = expect_pop(region);
        const double next_score = search_max_area(next_tree->copy_fixed_tree().get(), max_population) * pow(0.996, query_i + 1);
        if (query_i > 22 && current_score * pow(0.996, 1) > next_score)
            break;

        RegionNode* node_to_update = region_tree->find_equal_region_node(region);
        assert(node_to_update);
        assert(node_to_update->region == region);
        assert(!node_to_update->fixed);
        node_to_update->region.pop = query_region(region);
        node_to_update->fixed = true;
    }
//     return mark_predict_cities(world, prepre);

    unique_ptr<RegionNode> fixed_tree = region_tree->copy_fixed_tree();
    BitBoard selected = search_best_area_selection(fixed_tree.get(), max_population);
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
        g_timer.start();

//         return simple_search(max_percentage, world_map, total_population);

        const ll max_population = (ll)max_percentage * total_population / 100;
        BitBoard world(world_map);
        const ll total_area = world.count(0, 0, world.width(), world.height());

        BitBoard selected = select_area(max_population, world, total_population);
        fprintf(stderr, "run time: %f\n\n", g_timer.get_elapsed());
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
#ifdef CORRECT_POP
                    int p = correct_pop_map.pop(r);
#else
                    int p = query_region(x, y, r.high_x, r.high_y);
#endif
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
        vector<Rect> query_rect;
        if (smalls.size() - use_r.size() < use_r.size())
        {
            rep(i, smalls.size())
                if (!used_i[i])
                    query_rect.push_back(smalls[i]);
        }
        else
        {
            rep(i, smalls.size())
                if (used_i[i])
                    query_rect.push_back(smalls[i]);
        }
        dump(query_rect.size());

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

#ifdef CORRECT_POP
        for (auto& r : query_rect)
        {
            query_region(r.low_x, r.low_y, r.high_x, r.high_y);
        }
#endif

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
        return survey_result.selected;

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
        {
#ifdef CORRECT_POP
            rects[i].pop = correct_pop_map.pop(rects[i]);
#else
            rects[i].pop = query_region(rects[i].low_x, rects[i].low_y, rects[i].high_x, rects[i].high_y);
#endif
        }

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

//         fprintf(stderr, "my / opt: %8d %8d, %2.3f\n", best_area, opt_area, (double)best_area / opt_area);

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
    cout << "pop" << endl;
    cout.flush();

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

