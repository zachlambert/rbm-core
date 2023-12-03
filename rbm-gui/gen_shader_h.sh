# Usage:
# gen_shader_h.sh [shader file] [shader name] [parent include path]

dir="shaders/include/$3/shader"
mkdir -p $dir
path="$dir/$2.h"
echo "#ifndef GEN_SHADER_$2_h" > $path
echo "#define GEN_SHADER_$2_h" >> $path
echo "#include <string>" >> $path
echo "" >> $path
echo "namespace shader {" >> $path
echo "" >> $path
echo "static const std::string $2 = R\"V0G0N(" >> $path
cat shaders/$1 >> $path
echo ")V0G0N\";" >> $path
echo "" >> $path
echo "} //namespace shader" >> $path
echo "" >> $path
echo "#endif" >> $path
