(TeX-add-style-hook
 "ising_doc"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("scrartcl" "fleqn")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("natbib" "numbers")))
   (TeX-run-style-hooks
    "latex2e"
    "scrartcl"
    "scrartcl10"
    "amsmath"
    "amssymb"
    "bm"
    "graphics"
    "graphicx"
    "natbib")
   (TeX-add-symbols
    "D"
    "eul"
    "im"
    "tr")
   (LaTeX-add-labels
    "eq:can"
    "eq:balance"
    "eq:detbalance")
   (LaTeX-add-bibitems
    "newman99"
    "metropolis53")
   (LaTeX-add-environments
    "denseitem")))

