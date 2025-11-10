# Substitute Bars

Substitute the '+' function for the '\|' and '\|\|' function in a
mixed-model formula.

## Usage

``` r
subbars(term)
```

## Arguments

- term:

  a mixed-model formula

## Value

the formula with all \| and \|\| operators replaced by +

## Note

This function is called recursively on individual terms in the model,
which is why the argument is called `term` and not a name like `form`,
indicating a formula.
