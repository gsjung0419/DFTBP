# Install/uninstall DFTBP package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

mode=$1

# enforce using portable C locale
LC_ALL=C
export LC_ALL

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    cp $1 ..
    if (test $mode = 2) then
      echo "  updating src/$1"
    fi
  fi
}

for file in *.cpp *.h; do
  test -f "${file}" && action "${file}"
done

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*dftbp[^ \t]* //g' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(dftbp_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(dftbp_SYSLIB) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(dftbp_SYSPATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^[ \t]*include.*dftbp.*$/d' ../Makefile.package.settings
    sed -i -e '4 i \
include ..\/..\/lib\/dftbp\/Makefile.lammps
' ../Makefile.package.settings
  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^[ \t]*include.*dftbp.*$/d' ../Makefile.package.settings
  fi

fi
