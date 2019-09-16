build="build_sphinx/"
source="source/"
sphinx-build -b html $source $build
make html
echo "html file is now in $build"
