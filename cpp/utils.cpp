#include "utils.h"
#include <filesystem>

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
    fluid_dir = std::wstring(path.parent_path() / rel_fluid_dir);
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