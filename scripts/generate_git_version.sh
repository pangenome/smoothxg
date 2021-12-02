INC_DIR=$1

GIT_VERSION=$(git describe --always --tags)

echo "#define SMOOTHXG_GIT_VERSION" \"$GIT_VERSION\" > $1/smoothxg_git_version.hpp.tmp
diff $1/smoothxg_git_version.hpp.tmp $1/smoothxg_git_version.hpp >/dev/null || cp $1/smoothxg_git_version.hpp.tmp $1/smoothxg_git_version.hpp
rm -f $1/smoothxg_git_version.hpp.tmp
