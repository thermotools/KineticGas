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

    Value compute_if_absent(const Key& key, const std::function<Value(void)>& fun) {
        if (auto val = get(key)){
            return *val;
        }
        Value val = fun();
        return store_if_absent(key, val);
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

inline size_t combined_hash(const std::vector<size_t>& hashes) noexcept {
    // Pretty standard way to combine hash keys. See "hashing" and "golden ratio".
    size_t hash = 0;
    for (size_t h : hashes) {
        hash ^= h + 0x9e3779b9 + (hash << 6) + (hash >> 2);
    }
    return hash;
}

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
        return combined_hash({
            std::hash<int>{}(sp.T_dK),
            std::hash<double>{}(sp.rho)
        });
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
        int i = (op.i < op.j) ? op.i : op.j;
        int j = (op.i < op.j) ? op.j : op.i;
        return combined_hash({
            std::hash<int>{}(i),
            std::hash<int>{}(j),
            std::hash<int>{}(op.l),
            std::hash<int>{}(op.r),
            std::hash<int>{}(op.T_dK),
            std::hash<double>{}(op.rho)
        });
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
        int i = (s.i < s.j) ? s.i : s.j;
        int j = (s.i < s.j) ? s.j : s.i;
        return combined_hash({
            std::hash<int>{}(i),
            std::hash<int>{}(j),
            std::hash<int>{}(s.l),
            std::hash<double>{}(s.E)
        });
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