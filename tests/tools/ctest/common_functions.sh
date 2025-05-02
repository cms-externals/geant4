# common functions 

errexit()
{
  echo $1
  exit 1
}

warn()
{
  echo "Warning ================================================="
  echo "Warning : $1" 
  echo "Warning ================================================="
}

function testdir
{
  _name=$1
  _value=${!_name}
  test -z "$_value" && errexit "Error: variable $_name is not set" || true
  test -d $_value   && true || errexit "Error: directory for $_name = $_value not found"
}   

function PrependPath
{
  local _path=$1
  local _dir=$(echo $2 | sed -e 's/\/$//')  # drop trailing /
  if test -n "$_dir"; then 
    if echo ${!_path} | grep -q $_dir ; then
      true
    else
      test -n "${!_path}" \
        && export $_path=${_dir}:${!_path} \
        || export $_path=${_dir}
    fi      
  fi         
}

function AppendPath
{
  local _path=$1
  local _dir=$2
  if test -n "$_dir"; then 
    if echo ${!_path} | grep -q $_dir ; then
      true
    else
      test -n "${!_path}" \
        && export $_path=${!_path}:${_dir} \
        || export $_path=${_dir}
    fi
  fi         
}

function DeletePath
{
  local _path=$1
  local _dir=$(echo $2 | sed -e's/\/$//')
  if echo ${!_path} | grep -q $_dir ; then
    # break PATH into its elements and rebuild path from elements dropping _dir
    IFS=':' read -a elem <<<$(echo ${!_path})
    local newpath=''
    for ((i=0; i< ${#elem[*]};i++)) ; do
      if test "${elem[$i]}" = "${_dir}" -o "${elem[$i]}" = "${_dir}/" ; then
        # skip 
        true
      else
        # AppendPath would work, but also export newpath
        test -n "${newpath}" \
          && newpath=${newpath}:${elem[$i]} \
          || newpath=${elem[$i]}
      fi
    done
    export $_path=$newpath
  fi 
}

IsNotInPath() {
  package=$1
  rc=`echo ${PATH} | grep  $package | wc -l`
  return $rc
}

reduce_MAX_CPUS_USE () {
  test -z "$1" && return
  
  # value must be >= 1
  local result=$(( ( $1 >= 1 ) ? $1 : 1 ))
  
  # new value shall not be larger than current MAX_CPUS_USE, if defined
  test -n "$MAX_CPUS_USE" \
    && echo $(( ( $result <= $MAX_CPUS_USE ) ? $result : $MAX_CPUS_USE )) \
    ||  echo $result
}

