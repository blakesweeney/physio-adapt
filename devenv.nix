{
  pkgs,
  lib,
  config,
  inputs,
  ...
}:
{
  # https://devenv.sh/packages/
  packages = with pkgs; [
    claude-code
    easel
    git
    infernal
    nextflow
    ruff
    sqlite-utils
  ];

  # https://devenv.sh/languages/
  languages.python = {
    enable = true;
    venv.enable = true;
    uv.enable = true;
    uv.sync.enable = true;
  };

  # https://devenv.sh/git-hooks/
  git-hooks.hooks = {
    ruff.enable = true;
    nixfmt-rfc-style.enable = true;
  };

  # See full reference at https://devenv.sh/reference/options/
}
