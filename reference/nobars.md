# Omit terms separated by vertical bars in a formula

Remove the random-effects terms from a mixed-effects formula, thereby
producing the fixed-effects formula.

## Usage

``` r
nobars(term)
```

## Arguments

- term:

  the right-hand side of a mixed-model formula

## Value

the fixed-effects part of the formula

## Note

This function is called recursively on individual terms in the model,
which is why the argument is called `term` and not a name like `form`,
indicating a formula.
