#include <shared_mutex>
#include <optional>
#include <unordered_map>
#include <vector>

using vector2d = std::vector<std::vector<double>>;

template<typename Key, typename Value, typename KeyHash = std::hash<Key>>
class ThreadSafeMap {
public:
    bool has(const Key& key) const {
        std::shared_lock lock(mtx);
        return map.find(key) != map.end();
    }

    std::optional<Value> get(const Key& key) const {
        std::shared_lock lock(mtx);
        auto it = map.find(key);
        if (it != map.end()) {
            return it->second;
        }
        return std::nullopt;
    }

    void store(const Key& key, const Value& val){
        std::unique_lock lock(mtx);
        map[key] = val;
    }

    Value& store_if_absent(const Key& key, Value& val){
        {
            std::shared_lock lock(mtx);
            if (map.find(key) != map.end()){
                return val;
            }
        }
        std::unique_lock lock(mtx);
        if (map.find(key) != map.end()){
            return val;
        }
        map[key] = val;
        return val;
    }

    void clear(){
        std::unique_lock lock(mtx);
        map.clear();
    }
    
private:
    std::unordered_map<Key, Value, KeyHash> map;
    mutable std::shared_mutex mtx;
};

struct StatePoint{
    int T_dK;
    double rho;
    StatePoint(double T) : T_dK{static_cast<int>((T * 100.) + 0.5)}, rho{-1} {}
    StatePoint(double T, double rho) : T_dK{static_cast<int>((T * 100.) + 0.5)}, rho{rho}{}

    bool operator==(const StatePoint& other) const {
        return T_dK == other.T_dK && rho == other.rho;
    }
};

struct StatePointHash {
    size_t operator()(const StatePoint& sp) const noexcept {
        size_t h1 = std::hash<int>{}(sp.T_dK);
        size_t h2 = std::hash<double>{}(sp.rho);
        return h1 ^ (h2 << 1);
    }
};

struct OmegaPoint{
    int i, j, l, r, T_dK;
    double rho;
    OmegaPoint(int i, int j, int l, int r, double T, double rho) : i{i}, j{j}, l{l}, r{r}, rho{rho} {
         T_dK = (int) ((T * 10.0) + 0.5);
    };

    OmegaPoint(int i, int j, int l, int r, double T) : OmegaPoint(i, j, l, r, T, 0){}

    bool operator==(const OmegaPoint& other) const {
        return ((i == other.i && j == other.j) || (i == other.j && j == other.i))
                && l == other.l 
                && r == other.r
                && T_dK == other.T_dK 
                && rho == other.rho;
    }
};

struct OmegaPointHash {
    size_t operator()(const OmegaPoint& op) const noexcept {
        size_t h1 = std::hash<int>{}(op.i);
        size_t h2 = std::hash<int>{}(op.j);
        size_t h3 = std::hash<int>{}(op.l);
        size_t h4 = std::hash<int>{}(op.r);
        size_t h5 = std::hash<int>{}(op.T_dK);
        size_t h6 = std::hash<double>{}(op.rho);

        size_t seed = h1;
        seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h4 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h5 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= h6 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

struct CrossSectionPoint{
    int i, j, l, E;
    CrossSectionPoint(int i, int j, int l, double E)
        : i{i}, j{j}, l{l}, E{static_cast<int>(std::round(100 * E))} 
    {}

    bool operator==(const CrossSectionPoint& other) const{
        return ((i == other.i && j == other.j) || (i == other.j && j == other.i))
                && (l == other.l) 
                && (E == other.E);
    }
};

struct CrossSectionHash {
    size_t operator()(const CrossSectionPoint& s) const {
        return std::hash<int>()(s.i) ^ (std::hash<int>()(s.j) << 1) ^
               (std::hash<int>()(s.l) << 2) ^ (std::hash<double>()(s.E) << 3);
    }
};

class KineticGasCache {
public:
    void clear(){
        omega.clear();
        cross_section.clear();
        mtl.clear();
        etl.clear();
    }

    ThreadSafeMap<OmegaPoint, double, OmegaPointHash> omega;
    ThreadSafeMap<CrossSectionPoint, double, CrossSectionHash> cross_section;
    ThreadSafeMap<StatePoint, vector2d, StatePointHash> mtl;
    ThreadSafeMap<StatePoint, vector2d, StatePointHash> etl;
};