# Determine random-effects expressions from a formula

Determine random-effects expressions from a formula

## Usage

``` r
findbars(term)
```

## Arguments

- term:

  a mixed-model formula

## Value

pairs of expressions that were separated by vertical bars

## Note

This function is called recursively on individual terms in the model,
which is why the argument is called `term` and not a name like `form`,
indicating a formula.
