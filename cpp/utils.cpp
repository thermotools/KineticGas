#include "utils.h"
#include <filesystem>
#include <mutex>

#ifdef FLUID_DIR
    #define __MACRO_STRING__(s) __EXPAND_MACRO_STR__(s)
    #define __EXPAND_MACRO_STR__(s) #s
    static std::string rel_fluid_dir = __MACRO_STRING__(FLUID_DIR);
    #undef __MACRO_STRING__
    #undef __EXPAND_MACRO_STR__
#else
    static std::string rel_fluid_dir = "./fluids";
#endif

static std::string fluid_dir;

std::string get_fluid_dir(){
    if (fluid_dir.empty()) set_fluid_dir(rel_fluid_dir); // Ensure initialisation on first call
    return fluid_dir;
}

static void __set_fluid_dir(); // Platform specific handling of relative paths
void set_fluid_dir(const std::string path){
    rel_fluid_dir = path;
    __set_fluid_dir();
}

bool KineticGasCache::has_omega(const OmegaPoint& point) const {
    std::shared_lock lock(omega_mutex);
    return omega_map.find(point) != omega_map.end();
}

bool KineticGasCache::has_mtl(const StatePoint& stp) const {
    std::shared_lock lock(mtl_mutex);
    return mtl_map.find(stp) != mtl_map.end();
}

bool KineticGasCache::has_etl(const StatePoint& stp) const {
    std::shared_lock lock(etl_mutex);
    return etl_map.find(stp) != etl_map.end();
}

std::optional<double> KineticGasCache::get_omega(const OmegaPoint& point) const {
    std::shared_lock lock(omega_mutex);
    auto it = omega_map.find(point);
    if (it != omega_map.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::optional<vector2d> KineticGasCache::get_mtl(const StatePoint& stp) const {
    std::shared_lock lock(mtl_mutex);
    auto it = mtl_map.find(stp);
    if (it != mtl_map.end()) {
        return it->second;
    }
    return std::nullopt;
}

std::optional<vector2d> KineticGasCache::get_etl(const StatePoint& stp) const {
    std::shared_lock lock(etl_mutex);
    auto it = etl_map.find(stp);
    if (it != etl_map.end()) {
        return it->second;
    }
    return std::nullopt;
}

void KineticGasCache::store_omega(const OmegaPoint& point, const double omega) {
    std::unique_lock lock(omega_mutex);
    omega_map[point] = omega;
}

void KineticGasCache::store_mtl(const StatePoint& stp, const vector2d& mtl) {
    std::unique_lock lock(mtl_mutex);
    mtl_map[stp] = mtl;
}

void KineticGasCache::store_etl(const StatePoint& stp, const vector2d& etl) {
    std::unique_lock lock(etl_mutex);
    etl_map[stp] = etl;
}

void KineticGasCache::clear(){
    std::scoped_lock lock(omega_mutex, mtl_mutex, etl_mutex);
    omega_map.clear();
    mtl_map.clear();
    etl_map.clear();
}

#ifdef _WIN32
#include <windows.h>

static void __set_fluid_dir() {
    if (std::filesystem::path(rel_fluid_dir).is_absolute()){
        fluid_dir = rel_fluid_dir;
    }
    else {
    char chr_path[MAX_PATH];
    HMODULE hModule = NULL;
    GetModuleHandleEx(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT, (LPCSTR)&get_fluid_dir, &hModule);
    GetModuleFileName(hModule, chr_path, sizeof(chr_path));

    std::string libPath = std::string(chr_path);
    std::filesystem::path path(libPath);
    fluid_dir = (path.parent_path() / rel_fluid_dir).string();
    }
}

#elif defined(__APPLE__) || defined(__linux__)
#include <dlfcn.h>

static void __set_fluid_dir(){
    if (std::filesystem::path(rel_fluid_dir).is_absolute()){
        fluid_dir = rel_fluid_dir;
    }
    else {
        Dl_info info;
        dladdr((void*)&get_fluid_dir, &info);
        std::string libPath = info.dli_fname;
        std::filesystem::path path(libPath);
        fluid_dir = std::string(path.parent_path() / rel_fluid_dir);
    }
}

#else
    #error "Unknown platform!"
#endif