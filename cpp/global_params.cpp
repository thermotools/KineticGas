#include "global_params.h"

#define __MACRO_STRING__(s) __EXPAND_MACRO_STR__(s)
#define __EXPAND_MACRO_STR__(s) #s
#ifdef FLUID_DIR
    std::string fluid_dir = __MACRO_STRING__(FLUID_DIR);
#else
    std::string fluid_dir = "./fluids";
#endif
void set_fluid_dir(const std::string& path){
    fluid_dir = path;
}