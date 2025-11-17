# LT-FH++ Graph Example

``` r
library(LTFHPlus)
library(dplyr)
library(igraph)
```

### Automatic family construction

Manually constructing family history from registers can be tedious. It
can be easily done for first and second degree relatives, however, it
starts to become very difficult once further degrees of relatives needs
to be identified. Other R packages, such as *Kinship2* or *FamAgg*,
allows users to manipulate pedigrees and perform certain actions and
calculations on them, however, to the best of our knowledge, they do not
allow one to find all relatives of degree $`n`$ or closer of a given
individual. This particular problem is core to LT-FH++, as the family
members of a proband (often the genotyped individual) will need to be
identified as well as the family members’ mutual kinship coefficient for
the covariance matrix. With the introduction of the function
[`prepare_graph()`](https://emilmip.github.io/LTFHPlus/reference/prepare_graph.md),
we have made it possible to automatically construct families for a list
of probands. Below, we will provide an example of how the function can
be used in connection with the
[`estimate_liability()`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability.md)
function of the package. *All data is simulated and should not be used
as-is in any real-world analysis.*

### Constructing the population pedigree

The input is inspired by register information such as the trio
information provided by the danish CPR register. It usually follows the
format of having a column with the child’s ID (here simply `id`),
followed by columns for the mother’s and father’s ID (here `momcol` and
`dadcol`). Here, we provide an illustrative example of how the register
information may look. The names used are representative of their
relation to a proband named “pid”, with “pgm” refering to paternal
grandmother, “mgf” refering to maternal grandfather, the “hs” prefix
means half-sibling, the new parents of half-siblings will appear with
the “H” or “W” suffix for husband or wife, cousins are prefaced with the
first two letters of their parents, etc..

Additionally, we will also need threshold information for the
phenotype(s) of interest. The information is baked into the graph object
as attributes, and retrieved as needed later by the
[`estimate_liability()`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability.md)
function. It is possible to construct the graph without this information
and only construct the desired families.

``` r
# example of register trio input
family = tribble(
  ~id, ~momcol, ~dadcol,
  "pid", "mom", "dad",
  "sib", "mom", "dad",
  "mhs", "mom", "dad2",
  "phs", "mom2", "dad",
  "mom", "mgm", "mgf",
  "dad", "pgm", "pgf",
  "dad2", "pgm2", "pgf2",
  "paunt", "pgm", "pgf",
  "pacousin", "paunt", "pauntH",
  "hspaunt", "pgm", "newpgf",
  "hspacousin", "hspaunt", "hspauntH",
  "puncle", "pgm", "pgf",
  "pucousin", "puncleW", "puncle",
  "maunt", "mgm", "mgf",
  "macousin", "maunt", "mauntH",
  "hsmuncle", "newmgm", "mgf",
  "hsmucousin", "hsmuncleW", "hsmuncle"
) %>% print()
```

    ## # A tibble: 17 × 3
    ##    id         momcol    dadcol  
    ##    <chr>      <chr>     <chr>   
    ##  1 pid        mom       dad     
    ##  2 sib        mom       dad     
    ##  3 mhs        mom       dad2    
    ##  4 phs        mom2      dad     
    ##  5 mom        mgm       mgf     
    ##  6 dad        pgm       pgf     
    ##  7 dad2       pgm2      pgf2    
    ##  8 paunt      pgm       pgf     
    ##  9 pacousin   paunt     pauntH  
    ## 10 hspaunt    pgm       newpgf  
    ## 11 hspacousin hspaunt   hspauntH
    ## 12 puncle     pgm       pgf     
    ## 13 pucousin   puncleW   puncle  
    ## 14 maunt      mgm       mgf     
    ## 15 macousin   maunt     mauntH  
    ## 16 hsmuncle   newmgm    mgf     
    ## 17 hsmucousin hsmuncleW hsmuncle

``` r
# random thresholds to store as attributes
thrs =  tibble(
 id = family %>% select(1:3) %>% unlist() %>% unique(),
 lower = sample(c(-Inf, 2), size = length(id), replace = TRUE),
 upper = sample(c(2, Inf), size = length(id), replace = TRUE)) %>% print()
```

    ## # A tibble: 31 × 3
    ##    id       lower upper
    ##    <chr>    <dbl> <dbl>
    ##  1 pid       -Inf     2
    ##  2 sib       -Inf   Inf
    ##  3 mhs          2   Inf
    ##  4 phs       -Inf   Inf
    ##  5 mom       -Inf     2
    ##  6 dad          2     2
    ##  7 dad2      -Inf   Inf
    ##  8 paunt     -Inf     2
    ##  9 pacousin  -Inf     2
    ## 10 hspaunt      2   Inf
    ## # ℹ 21 more rows

### Constructing the graph

With the above input, we can construct the graph for the population
specified in `family` and attach the information in `thrs`:

``` r
graph = prepare_graph(.tbl = family, 
                      node_attributes = thrs,
                      fcol = "dadcol",
                      mcol = "momcol",
                      icol = "id")
graph
```

    ## IGRAPH 750e52d DN-- 31 44 -- 
    ## + attr: name (v/c), lower (v/n), upper (v/n)
    ## + edges from 750e52d (vertex names):
    ##  [1] dad     ->pid        mom     ->pid        dad     ->sib       
    ##  [4] mom     ->sib        dad2    ->mhs        mom     ->mhs       
    ##  [7] dad     ->phs        mom2    ->phs        mgf     ->mom       
    ## [10] mgm     ->mom        pgf     ->dad        pgm     ->dad       
    ## [13] pgf2    ->dad2       pgm2    ->dad2       pgf     ->paunt     
    ## [16] pgm     ->paunt      pauntH  ->pacousin   paunt   ->pacousin  
    ## [19] newpgf  ->hspaunt    pgm     ->hspaunt    hspauntH->hspacousin
    ## [22] hspaunt ->hspacousin pgf     ->puncle     pgm     ->puncle    
    ## + ... omitted several edges

The object `graph` is a directed graph with the connections specified in
the `family` object and attributes (attr) from `thrs`. The package
*igraph* can be used to do any manipulation of the object `graph`. The
object `graph` only contains a single family in this example, however,
it could just as easily contain the entirety of the danish CPR register,
totalling more than $`10`$ million individuals. We can create a local
graph (formally, a neighborhood sub-graph) around a proband, only
including individuals that are within $`n`$ steps of said proband. The
$`n`$ steps correspond to $`n`$ degree relatives, such that $`n = 1`$
yields all first degree relatives, $`n = 2`$ yields all second degree
relatives, etc.. From this local graph, we can construct a kinship
matrix and the covariance matrix required by
[`estimate_liability()`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability.md).

### Extracting local graphs

If we have constructed a graph on the entire Danish CPR register and a
subset of these individuals have been genotyped, then we can create
local graphs around each of the genotyped individuals. Here, we simply
use every individual for the illustration:

``` r
genotyped_ids = V(graph)$name

family_graphs = tibble(
  pid = genotyped_ids,
  fam_graph = make_ego_graph(graph, order = 2, nodes = pid)
)
```

Here, `family_graphs` is a tibble with an id column called “pid” and a
column called “fam_graph” with list of all the local graphs, centred on
the individual from “pid”.

### Construct kinship matrix

We can construct a kinship matrix for a local graph in the following
way:

``` r
get_kinship(family_graphs$fam_graph[[1]], h2 = 1, index_id = family_graphs$pid[1])
```

    ##         pid  sib  mhs  phs mom dad paunt puncle maunt  mgm  pgm  mgf  pgf pid_g
    ## pid    1.00 0.50 0.25 0.25 0.5 0.5  0.25   0.25  0.25 0.25 0.25 0.25 0.25  1.00
    ## sib    0.50 1.00 0.25 0.25 0.5 0.5  0.25   0.25  0.25 0.25 0.25 0.25 0.25  0.50
    ## mhs    0.25 0.25 1.00 0.00 0.5 0.0  0.00   0.00  0.25 0.25 0.00 0.25 0.00  0.25
    ## phs    0.25 0.25 0.00 1.00 0.0 0.5  0.25   0.25  0.00 0.00 0.25 0.00 0.25  0.25
    ## mom    0.50 0.50 0.50 0.00 1.0 0.0  0.00   0.00  0.50 0.50 0.00 0.50 0.00  0.50
    ## dad    0.50 0.50 0.00 0.50 0.0 1.0  0.50   0.50  0.00 0.00 0.50 0.00 0.50  0.50
    ## paunt  0.25 0.25 0.00 0.25 0.0 0.5  1.00   0.50  0.00 0.00 0.50 0.00 0.50  0.25
    ## puncle 0.25 0.25 0.00 0.25 0.0 0.5  0.50   1.00  0.00 0.00 0.50 0.00 0.50  0.25
    ## maunt  0.25 0.25 0.25 0.00 0.5 0.0  0.00   0.00  1.00 0.50 0.00 0.50 0.00  0.25
    ## mgm    0.25 0.25 0.25 0.00 0.5 0.0  0.00   0.00  0.50 1.00 0.00 0.00 0.00  0.25
    ## pgm    0.25 0.25 0.00 0.25 0.0 0.5  0.50   0.50  0.00 0.00 1.00 0.00 0.00  0.25
    ## mgf    0.25 0.25 0.25 0.00 0.5 0.0  0.00   0.00  0.50 0.00 0.00 1.00 0.00  0.25
    ## pgf    0.25 0.25 0.00 0.25 0.0 0.5  0.50   0.50  0.00 0.00 0.00 0.00 1.00  0.25
    ## pid_g  1.00 0.50 0.25 0.25 0.5 0.5  0.25   0.25  0.25 0.25 0.25 0.25 0.25  1.00

Note: individuals such as “hspaunt” are not present, since they are
third degree relation.

### Calculating genetic liabilities

With the object `family_graphs`, we have constructed an object with the
required information and format for estimating the genetic liabilities
of a list of probands. We have automatically extracted the family
present for a given proband up to degree $`2`$ (can be generalised to
any integer $`n`$), which allows us to construct the kinship matrix.
Next, we have attached the lower and upper threshold for each individual
during the construction of the graph, and as such follow the indivdual
in whatever local graph they are a part of. Finally, we need a
liability-scale heritability. We will simply use $`h^2 = 0.5`$ for
illustrative purposes. We will use
[`estimate_liability()`](https://emilmip.github.io/LTFHPlus/reference/estimate_liability.md)
in the folling way:

``` r
ltfhpp = estimate_liability(
  family_graphs = family_graphs,
  h2 = .5,
  pid = "pid",
  family_graphs_col = "fam_graph"
) %>% print()
```

    ## The number of workers is 1

    ## # A tibble: 31 × 4
    ##    fam_ID   pid      genetic_est genetic_se
    ##    <chr>    <chr>          <dbl>      <dbl>
    ##  1 pid      pid            0.939    0.00270
    ##  2 sib      sib            1.01     0.00350
    ##  3 mhs      mhs            1.44     0.00208
    ##  4 phs      phs            0.566    0.00349
    ##  5 mom      mom            1.07     0.00284
    ##  6 dad      dad            1.21     0.00186
    ##  7 dad2     dad2           1.07     0.00341
    ##  8 paunt    paunt          0.773    0.00272
    ##  9 pacousin pacousin       0.301    0.00366
    ## 10 hspaunt  hspaunt        1.85     0.00193
    ## # ℹ 21 more rows
