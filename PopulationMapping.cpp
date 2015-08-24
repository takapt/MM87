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
        T v = next_double(sum) + (T)1e-6;
        for (int i = 0; i < (int)ratio.size(); ++i)
        {
            v -= ratio[i];
            if (v <= 0)
                return i;
        }
        return 0;
    }
};
Random g_rand;


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

class World
{
public:
    World(const vector<string>& s) :
        h(s.size()), w(s[0].size())
    {
        rep(y, h) rep(x, w)
            if (s[y][x] == 'X')
                f[index(x, y)] = true;
    }

    bool is_land(int x, int y) const
    {
        return f[index(x, y)];
    }

    int area(int low_x, int low_y, int high_x, int high_y) const
    {
        int c = 0;
        for (int y = low_y; y < high_y; ++y)
            for (int x = low_x; x < high_x; ++x)
                if (is_land(x, y))
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

    void move(int dx, int dy)
    {
        low_x += dx;
        high_x += dx;
        low_y += dy;
        high_y += dy;
    }

    // true if either region contains the other
    bool intersect(const Region& other) const
    {
        return low_x < other.high_x && other.low_x < high_x &&
               low_y < other.high_y && other.low_y < high_y;
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
        return copy_node;
    }

    const vector<unique_ptr<RegionNode>>& childs_nodes() const
    {
        return childs;
    }

    void add_child_region(const Region& reg)
    {
        assert(region.contain(reg));
        assert(!intersect_with_child(reg));
        childs.push_back(unique_ptr<RegionNode>(new RegionNode(reg, this, false)));
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
        assert(area >= 0);
        return area;
    }

    ll pop_except_childs() const
    {
        int pop = region.pop;
        for (auto& child : childs)
            pop -= child->region.pop;
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

    int tree_size() const
    {
        int size = 0;
        iterate([&](const RegionNode*) { ++size; });
        return size;
    }

    vector<RegionNode*> list_all(bool ignore_fixed)
    {
        vector<RegionNode*> nodes;
        iterate([&](RegionNode* no) {
            if (!fixed || !ignore_fixed)
                nodes.push_back(no);
        });
        return nodes;
    }
    vector<RegionNode*> list_leaves(bool ignore_fixed)
    {
        vector<RegionNode*> nodes;
        iterate([&](RegionNode* no) {
            if (no->childs.empty() && (!fixed || !ignore_fixed))
                nodes.push_back(no);
        });
        return nodes;
    }

    Region region;
    RegionNode* parent;

    bool fixed = false;

private:
    vector<unique_ptr<RegionNode>> childs;
};

struct City
{
    City(int x, int y, int scale) :
        x(x), y(y), scale(scale)
    {
        assert(1 <= scale && scale <= 5);
    }

    int x, y, scale;
};

class PopMap
{
public:
    PopMap(int w, int h) :
        w(w), h(h), p{}
    {
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
#define CORRECT_POP
#ifdef CORRECT_POP
PopMap correct_pop_map;
#endif

vector<ll> search_min_pop_for_area(const RegionNode* region_tree)
{
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

vector<Region> search_queries(const World& world, const PopMap& pop_map, const ll max_population, const RegionNode* fixed_region_tree)
{
    const auto generate_region = [&](const RegionNode* tree_) -> unique_ptr<RegionNode>
    {
        auto tree = tree_->copy_tree();
        auto nodes = tree->list_all(true);
        if (nodes.empty())
            return nullptr;

        RegionNode* target_region_node = nodes[g_rand.next_int(nodes.size())];
        const Region& target = target_region_node->region;

        int low_x = g_rand.next_int(target.low_x, target.high_x);
        int low_y = g_rand.next_int(target.low_y, target.high_y);
        int high_x = g_rand.next_int(low_x, target.high_x) + 1;
        int high_y = g_rand.next_int(low_y, target.high_y) + 1;
        Region added_region(low_x, low_y, high_x, high_y);
        assert(target.contain(added_region));
        if (target_region_node->intersect_with_child(added_region))
            return nullptr;

        added_region.area = world.area(low_x, low_y, high_x, high_y);
        added_region.pop = pop_map.pop(low_x, low_y, high_x, high_y);
        target_region_node->add_child_region(added_region);

        return tree;
    };

    const auto remove_region = [](const RegionNode* tree_) -> unique_ptr<RegionNode>
    {
        auto tree = tree_->copy_tree();
        auto nodes = tree->list_leaves(true);
        if (nodes.empty())
            return nullptr;

        RegionNode* remove_node = nodes[g_rand.next_int(nodes.size())];
        remove_node->parent->remove_child_region(remove_node);

        return tree;
    };

    const auto move_region = [](const RegionNode* tree_) -> unique_ptr<RegionNode>
    {
        auto tree = tree_->copy_tree();
        auto nodes = tree->list_all(true);
        if (nodes.empty())
            return nullptr;

        RegionNode* move_region_node = nodes[g_rand.next_int(nodes.size())];
        const RegionNode* parent = move_region_node->parent;
        const Region& parent_region = parent->region;

        Region move_region = move_region_node->region;
        // [low, high]
        int low_dx = -20, high_dx = 20, low_dy = -20, high_dy = 20;
        upmax(low_dx, parent_region.low_x - move_region.low_x);
        upmin(high_dx, parent_region.high_x - move_region.high_x);
        upmax(low_dy, parent_region.low_y - move_region.low_y);
        upmin(high_dy, parent_region.high_y - move_region.high_y);
        for (auto& child : move_region_node->childs_nodes())
        {
            const Region& reg = child->region;
            upmax(low_dx, move_region.high_x - reg.high_x);
            upmin(high_dx, move_region.low_x - reg.low_x);
            upmax(low_dy, move_region.high_y - reg.high_y);
            upmin(high_dy, move_region.low_y - reg.low_y);
        }
        assert(low_dx <= 0 && 0 <= high_dx);
        assert(low_dy <= 0 && 0 <= high_dy);
        assert(low_dx <= high_dx);
        assert(low_dy <= high_dy);

        int dx = g_rand.next_int(low_dx, high_dx + 1);
        int dy = g_rand.next_int(low_dy, high_dy + 1);
        move_region.move(dx, dy);
        assert(!move_region_node->parent->region.intersect_on_side(move_region));
        assert(!move_region_node->intersect_with_child(move_region));
        if (parent->intersect_with_child(move_region))
            return nullptr;

        move_region_node->region = move_region;
        return tree;
    };

    const auto resize_region = [](const RegionNode* tree_) -> unique_ptr<RegionNode>
    {
        auto tree = tree_->copy_tree();
        auto nodes = tree->list_all(true);
        if (nodes.empty())
            return nullptr;

        RegionNode* resize_region_node = nodes[g_rand.next_int(nodes.size())];
        const RegionNode* parent = resize_region_node->parent;
        const Region& parent_region = parent->region;

        Region resize_region = resize_region_node->region;
        // [low, high]
        int low_dx = max(-20, -(resize_region.high_x - resize_region.low_x - 1));
        int high_dx = min(20, parent_region.high_x - resize_region.high_x);
        int low_dy = max(-20, -(resize_region.high_y - resize_region.low_y - 1));
        int high_dy = min(20, parent_region.high_y - resize_region.high_y);
        for (auto& child : resize_region_node->childs_nodes())
        {
            const Region& reg = child->region;
            upmax(low_dx, resize_region.high_x - reg.high_x);
            upmax(low_dy, resize_region.high_y - reg.high_y);
        }
        assert(low_dx <= 0 && 0 <= high_dx);
        assert(low_dy <= 0 && 0 <= high_dy);
        assert(low_dx <= high_dx);
        assert(low_dy <= high_dy);

        int dx = g_rand.next_int(low_dx, high_dx + 1);
        int dy = g_rand.next_int(low_dy, high_dy + 1);
        resize_region.high_x += dx;
        resize_region.high_y += dy;
        assert(resize_region.low_x < resize_region.high_x);
        assert(resize_region.low_y < resize_region.high_y);
        assert(!parent_region.intersect_on_side(resize_region));
        assert(!resize_region_node->intersect_with_child_on_side(resize_region));
        if ((dx >= 0 || dy >= 0) && parent->intersect_with_child(resize_region))
            return nullptr;

        resize_region_node->region = resize_region;
        return tree;
    };

    unique_ptr<RegionNode> current_tree = fixed_region_tree->copy_tree();
    current_tree->iterate([](RegionNode* no) { no->fixed = true; });

    return {};
}


class PopulationMapping
{
public:
    vector<string> mapPopulation(int max_percentage, vector <string> world_map, int total_population)
    {
        const int max_population = (ll)total_population * max_percentage / 100;
        const int h = world_map.size(), w = world_map[0].size();
        vector<vector<bool>> selected(h, vector<bool>(w));

        World world(world_map);

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
                rect.area = world.area(lx, ly, hx, hy);
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

private:
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
        World world(world_map);

        vector<Rect> smalls;
        const int S = survey_size;
        for (int y = 0; y < h; y += S)
        {
            for (int x = 0; x < w; x += S)
            {
                Rect r(x, y, min(x + S, w), min(y + S, h));
                r.area = world.area(x, y, r.high_x, r.high_y);
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

        const int total_area = world.area(0, 0, world.width(), world.height());

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
#endif
};


#ifdef LOCAL
int main(int argc, char* argv[])
{
    if (argc > 1)
        survey_size = to_T<int>(argv[1]);

    int max_percentage;
    cin >> max_percentage;
    int height;
    cin >> height;
    vector<string> world_map(height);
    input(world_map, height);
    int total_population;
    cin >> total_population;

#ifdef CORRECT_POP
    int h, w;
    cin >> h >> w;
    correct_pop_map = PopMap(w, h);
    rep(y, h) rep(x, w)
    {
        int p;
        cin >> p;
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
