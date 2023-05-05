alias python="python3.9"
alias exit_container="exit"
alias so="source ~/.bashrc"

function build(){
    cd /root/code
    rm -Rf cpp/release/*
    cd cpp/release && cmake -DCMAKE_BUILD_TYPE=Release  .. && make
    cp cpp/release/KineticGas*.so pykingas/KineticGas.so
}

alias b="build"
