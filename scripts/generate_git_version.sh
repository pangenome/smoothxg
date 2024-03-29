INC_DIR=$1

# Go to the directory where the script is
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd "$SCRIPT_DIR"

GIT_VERSION=$(git describe --always --tags --long)

echo "#define SMOOTHXG_GIT_VERSION" \"$GIT_VERSION\" > "$INC_DIR"/smoothxg_git_version.hpp
