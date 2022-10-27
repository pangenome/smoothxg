INC_DIR=$1

GIT_VERSION=$(git describe --always --tags)

echo "#define SMOOTHXG_GIT_VERSION" \"$GIT_VERSION\" > "$INC_DIR"/smoothxg_git_version.hpp
