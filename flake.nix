{
  description = "Unfolding Muon misID samples.";

  inputs = {
    root-curated.url = "github:umd-lhcb/root-curated";
    nixpkgs.follows = "root-curated/nixpkgs";
    flake-utils.follows = "root-curated/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils, root-curated }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config = { allowUnfree = true; };
          overlays = [ root-curated.overlay ];
        };
        python = pkgs.python3;
        pythonPackages = python.pkgs;
      in
      rec {
        packages = flake-utils.lib.flattenTree {
          dev-shell = devShell.inputDerivation;
        };
        devShell = pkgs.mkShell rec {
          name = "misid-unfold-dev";
          buildInputs = (with pkgs; with pythonPackages; [
            # Dev tools
            clang-tools

            root
            cxxopts
            roounfold

            # Python stack
            #pythonPackages.pyTuplingUtils
            virtualenvwrapper
            pylint

            # Pinned Python dependencies
            numpy
          ]);

          FONTCONFIG_FILE = pkgs.makeFontsConf {
            fontDirectories = with pkgs; [
              corefonts
            ];
          };

          shellHook = ''
            export PATH=$(pwd)/bin:$PATH

            # Allow the use of wheels.
            SOURCE_DATE_EPOCH=$(date +%s)

            VENV=./.virtualenv

            if test ! -d $VENV; then
              virtualenv $VENV
            fi
            source $VENV/bin/activate

            # allow for the environment to pick up packages installed with virtualenv
            export PYTHONPATH=$VENV/${python.sitePackages}/:$PYTHONPATH

            # fix libstdc++.so not found error
            export LD_LIBRARY_PATH=${pkgs.stdenv.cc.cc.lib}/lib:$LD_LIBRARY_PATH
          '';
        };
      });
}
