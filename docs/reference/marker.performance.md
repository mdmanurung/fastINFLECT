# Plot marker performance across metaclustering results

Plots criterion-specific marker QC across metaclusterings. The score is
the percentage of clusters where a marker passed the aggregate criterion
selected in `INFLECT()`. It is a screening pass rate, not proof of
unimodality. The colour gradient denotes the tested cluster count.

## Usage

``` r
marker.performance(inflect.results, ggtitle = NULL, markers = NULL)
```

## Arguments

  - inflect.results:
    
    A inflect.results object resulting from `INFLECT`.

  - ggtitle:
    
    Optional. Character string to plot as title.

  - markers:
    
    Which markers should be included in the plot? A vector of strings
    matching evaluated marker names. If `NULL`, all evaluated markers
    are displayed.

## Value

A list with `marker.dataframe` and `plot`. The data frame contains
canonical `k`, `marker`, and `qc_pass_rate` columns. Deprecated `i`,
`Marker`, and `Performance` aliases are retained for compatibility.

## Examples

``` r
# Read in FlowSOM object from file. Downsampled clustering result of Levine32 dataset clustering.
# SOM-clustered to 375 clusters.
flowsom <- system.file("extdata", "Levine32sample.Rdata", package="fastINFLECT")
load(flowsom)
inflect.results <- INFLECT(
  FlowSOM.results = dataset,
  set.i = 5:12,
  multicore = FALSE
)

# Display diagnostic graph
inflect.results$ggplot


# Now check marker performance for all markers
marker.performance(inflect.results, ggtitle= "Levine32sample", markers=NULL)
#> $marker.dataframe
#>      k           marker qc_pass_rate  i           Marker Performance
#> 1    5     CD3(Er170)Di     20.00000  5     CD3(Er170)Di    20.00000
#> 2    6     CD3(Er170)Di     50.00000  6     CD3(Er170)Di    50.00000
#> 3    7     CD3(Er170)Di     57.14286  7     CD3(Er170)Di    57.14286
#> 4    8     CD3(Er170)Di     62.50000  8     CD3(Er170)Di    62.50000
#> 5    9     CD3(Er170)Di     55.55556  9     CD3(Er170)Di    55.55556
#> 6   10     CD3(Er170)Di     60.00000 10     CD3(Er170)Di    60.00000
#> 7   11     CD3(Er170)Di     63.63636 11     CD3(Er170)Di    63.63636
#> 8   12     CD3(Er170)Di     66.66667 12     CD3(Er170)Di    66.66667
#> 9    5     CD4(Nd145)Di     40.00000  5     CD4(Nd145)Di    40.00000
#> 10   6     CD4(Nd145)Di     50.00000  6     CD4(Nd145)Di    50.00000
#> 11   7     CD4(Nd145)Di     71.42857  7     CD4(Nd145)Di    71.42857
#> 12   8     CD4(Nd145)Di     62.50000  8     CD4(Nd145)Di    62.50000
#> 13   9     CD4(Nd145)Di     66.66667  9     CD4(Nd145)Di    66.66667
#> 14  10     CD4(Nd145)Di     70.00000 10     CD4(Nd145)Di    70.00000
#> 15  11     CD4(Nd145)Di     72.72727 11     CD4(Nd145)Di    72.72727
#> 16  12     CD4(Nd145)Di     75.00000 12     CD4(Nd145)Di    75.00000
#> 17   5     CD7(Dy162)Di     80.00000  5     CD7(Dy162)Di    80.00000
#> 18   6     CD7(Dy162)Di    100.00000  6     CD7(Dy162)Di   100.00000
#> 19   7     CD7(Dy162)Di    100.00000  7     CD7(Dy162)Di   100.00000
#> 20   8     CD7(Dy162)Di    100.00000  8     CD7(Dy162)Di   100.00000
#> 21   9     CD7(Dy162)Di    100.00000  9     CD7(Dy162)Di   100.00000
#> 22  10     CD7(Dy162)Di    100.00000 10     CD7(Dy162)Di   100.00000
#> 23  11     CD7(Dy162)Di    100.00000 11     CD7(Dy162)Di   100.00000
#> 24  12     CD7(Dy162)Di    100.00000 12     CD7(Dy162)Di   100.00000
#> 25   5     CD8(Nd146)Di     20.00000  5     CD8(Nd146)Di    20.00000
#> 26   6     CD8(Nd146)Di     33.33333  6     CD8(Nd146)Di    33.33333
#> 27   7     CD8(Nd146)Di     57.14286  7     CD8(Nd146)Di    57.14286
#> 28   8     CD8(Nd146)Di     62.50000  8     CD8(Nd146)Di    62.50000
#> 29   9     CD8(Nd146)Di     66.66667  9     CD8(Nd146)Di    66.66667
#> 30  10     CD8(Nd146)Di     70.00000 10     CD8(Nd146)Di    70.00000
#> 31  11     CD8(Nd146)Di     72.72727 11     CD8(Nd146)Di    72.72727
#> 32  12     CD8(Nd146)Di     75.00000 12     CD8(Nd146)Di    75.00000
#> 33   5   CD11b(Nd144)Di     20.00000  5   CD11b(Nd144)Di    20.00000
#> 34   6   CD11b(Nd144)Di     33.33333  6   CD11b(Nd144)Di    33.33333
#> 35   7   CD11b(Nd144)Di     28.57143  7   CD11b(Nd144)Di    28.57143
#> 36   8   CD11b(Nd144)Di     25.00000  8   CD11b(Nd144)Di    25.00000
#> 37   9   CD11b(Nd144)Di     33.33333  9   CD11b(Nd144)Di    33.33333
#> 38  10   CD11b(Nd144)Di     40.00000 10   CD11b(Nd144)Di    40.00000
#> 39  11   CD11b(Nd144)Di     45.45455 11   CD11b(Nd144)Di    45.45455
#> 40  12   CD11b(Nd144)Di     50.00000 12   CD11b(Nd144)Di    50.00000
#> 41   5   CD11c(Tb159)Di     80.00000  5   CD11c(Tb159)Di    80.00000
#> 42   6   CD11c(Tb159)Di     83.33333  6   CD11c(Tb159)Di    83.33333
#> 43   7   CD11c(Tb159)Di     85.71429  7   CD11c(Tb159)Di    85.71429
#> 44   8   CD11c(Tb159)Di     87.50000  8   CD11c(Tb159)Di    87.50000
#> 45   9   CD11c(Tb159)Di     88.88889  9   CD11c(Tb159)Di    88.88889
#> 46  10   CD11c(Tb159)Di     90.00000 10   CD11c(Tb159)Di    90.00000
#> 47  11   CD11c(Tb159)Di     90.90909 11   CD11c(Tb159)Di    90.90909
#> 48  12   CD11c(Tb159)Di     91.66667 12   CD11c(Tb159)Di    91.66667
#> 49   5    CD13(Er168)Di     80.00000  5    CD13(Er168)Di    80.00000
#> 50   6    CD13(Er168)Di     83.33333  6    CD13(Er168)Di    83.33333
#> 51   7    CD13(Er168)Di     85.71429  7    CD13(Er168)Di    85.71429
#> 52   8    CD13(Er168)Di     87.50000  8    CD13(Er168)Di    87.50000
#> 53   9    CD13(Er168)Di     88.88889  9    CD13(Er168)Di    88.88889
#> 54  10    CD13(Er168)Di     90.00000 10    CD13(Er168)Di    90.00000
#> 55  11    CD13(Er168)Di     90.90909 11    CD13(Er168)Di    90.90909
#> 56  12    CD13(Er168)Di     91.66667 12    CD13(Er168)Di    91.66667
#> 57   5    CD14(Gd156)Di    100.00000  5    CD14(Gd156)Di   100.00000
#> 58   6    CD14(Gd156)Di    100.00000  6    CD14(Gd156)Di   100.00000
#> 59   7    CD14(Gd156)Di    100.00000  7    CD14(Gd156)Di   100.00000
#> 60   8    CD14(Gd156)Di    100.00000  8    CD14(Gd156)Di   100.00000
#> 61   9    CD14(Gd156)Di    100.00000  9    CD14(Gd156)Di   100.00000
#> 62  10    CD14(Gd156)Di    100.00000 10    CD14(Gd156)Di   100.00000
#> 63  11    CD14(Gd156)Di    100.00000 11    CD14(Gd156)Di   100.00000
#> 64  12    CD14(Gd156)Di    100.00000 12    CD14(Gd156)Di   100.00000
#> 65   5    CD15(Dy164)Di     60.00000  5    CD15(Dy164)Di    60.00000
#> 66   6    CD15(Dy164)Di     66.66667  6    CD15(Dy164)Di    66.66667
#> 67   7    CD15(Dy164)Di     71.42857  7    CD15(Dy164)Di    71.42857
#> 68   8    CD15(Dy164)Di     87.50000  8    CD15(Dy164)Di    87.50000
#> 69   9    CD15(Dy164)Di     88.88889  9    CD15(Dy164)Di    88.88889
#> 70  10    CD15(Dy164)Di     90.00000 10    CD15(Dy164)Di    90.00000
#> 71  11    CD15(Dy164)Di     90.90909 11    CD15(Dy164)Di    90.90909
#> 72  12    CD15(Dy164)Di    100.00000 12    CD15(Dy164)Di   100.00000
#> 73   5    CD16(Ho165)Di     40.00000  5    CD16(Ho165)Di    40.00000
#> 74   6    CD16(Ho165)Di     50.00000  6    CD16(Ho165)Di    50.00000
#> 75   7    CD16(Ho165)Di     71.42857  7    CD16(Ho165)Di    71.42857
#> 76   8    CD16(Ho165)Di     62.50000  8    CD16(Ho165)Di    62.50000
#> 77   9    CD16(Ho165)Di     66.66667  9    CD16(Ho165)Di    66.66667
#> 78  10    CD16(Ho165)Di     70.00000 10    CD16(Ho165)Di    70.00000
#> 79  11    CD16(Ho165)Di     63.63636 11    CD16(Ho165)Di    63.63636
#> 80  12    CD16(Ho165)Di     66.66667 12    CD16(Ho165)Di    66.66667
#> 81   5    CD19(Nd142)Di     40.00000  5    CD19(Nd142)Di    40.00000
#> 82   6    CD19(Nd142)Di     50.00000  6    CD19(Nd142)Di    50.00000
#> 83   7    CD19(Nd142)Di     57.14286  7    CD19(Nd142)Di    57.14286
#> 84   8    CD19(Nd142)Di     50.00000  8    CD19(Nd142)Di    50.00000
#> 85   9    CD19(Nd142)Di     55.55556  9    CD19(Nd142)Di    55.55556
#> 86  10    CD19(Nd142)Di     60.00000 10    CD19(Nd142)Di    60.00000
#> 87  11    CD19(Nd142)Di     72.72727 11    CD19(Nd142)Di    72.72727
#> 88  12    CD19(Nd142)Di     75.00000 12    CD19(Nd142)Di    75.00000
#> 89   5    CD20(Sm147)Di     40.00000  5    CD20(Sm147)Di    40.00000
#> 90   6    CD20(Sm147)Di     50.00000  6    CD20(Sm147)Di    50.00000
#> 91   7    CD20(Sm147)Di     71.42857  7    CD20(Sm147)Di    71.42857
#> 92   8    CD20(Sm147)Di     62.50000  8    CD20(Sm147)Di    62.50000
#> 93   9    CD20(Sm147)Di     66.66667  9    CD20(Sm147)Di    66.66667
#> 94  10    CD20(Sm147)Di     70.00000 10    CD20(Sm147)Di    70.00000
#> 95  11    CD20(Sm147)Di     72.72727 11    CD20(Sm147)Di    72.72727
#> 96  12    CD20(Sm147)Di     75.00000 12    CD20(Sm147)Di    75.00000
#> 97   5    CD22(Nd143)Di     40.00000  5    CD22(Nd143)Di    40.00000
#> 98   6    CD22(Nd143)Di     50.00000  6    CD22(Nd143)Di    50.00000
#> 99   7    CD22(Nd143)Di     57.14286  7    CD22(Nd143)Di    57.14286
#> 100  8    CD22(Nd143)Di     50.00000  8    CD22(Nd143)Di    50.00000
#> 101  9    CD22(Nd143)Di     55.55556  9    CD22(Nd143)Di    55.55556
#> 102 10    CD22(Nd143)Di     60.00000 10    CD22(Nd143)Di    60.00000
#> 103 11    CD22(Nd143)Di     72.72727 11    CD22(Nd143)Di    72.72727
#> 104 12    CD22(Nd143)Di     75.00000 12    CD22(Nd143)Di    75.00000
#> 105  5    CD33(Gd158)Di     40.00000  5    CD33(Gd158)Di    40.00000
#> 106  6    CD33(Gd158)Di     50.00000  6    CD33(Gd158)Di    50.00000
#> 107  7    CD33(Gd158)Di     57.14286  7    CD33(Gd158)Di    57.14286
#> 108  8    CD33(Gd158)Di     50.00000  8    CD33(Gd158)Di    50.00000
#> 109  9    CD33(Gd158)Di     55.55556  9    CD33(Gd158)Di    55.55556
#> 110 10    CD33(Gd158)Di     60.00000 10    CD33(Gd158)Di    60.00000
#> 111 11    CD33(Gd158)Di     63.63636 11    CD33(Gd158)Di    63.63636
#> 112 12    CD33(Gd158)Di     66.66667 12    CD33(Gd158)Di    66.66667
#> 113  5    CD34(Nd148)Di     20.00000  5    CD34(Nd148)Di    20.00000
#> 114  6    CD34(Nd148)Di     16.66667  6    CD34(Nd148)Di    16.66667
#> 115  7    CD34(Nd148)Di     14.28571  7    CD34(Nd148)Di    14.28571
#> 116  8    CD34(Nd148)Di     12.50000  8    CD34(Nd148)Di    12.50000
#> 117  9    CD34(Nd148)Di     22.22222  9    CD34(Nd148)Di    22.22222
#> 118 10    CD34(Nd148)Di     30.00000 10    CD34(Nd148)Di    30.00000
#> 119 11    CD34(Nd148)Di     45.45455 11    CD34(Nd148)Di    45.45455
#> 120 12    CD34(Nd148)Di     50.00000 12    CD34(Nd148)Di    50.00000
#> 121  5    CD38(Er167)Di     60.00000  5    CD38(Er167)Di    60.00000
#> 122  6    CD38(Er167)Di     50.00000  6    CD38(Er167)Di    50.00000
#> 123  7    CD38(Er167)Di     71.42857  7    CD38(Er167)Di    71.42857
#> 124  8    CD38(Er167)Di     75.00000  8    CD38(Er167)Di    75.00000
#> 125  9    CD38(Er167)Di     77.77778  9    CD38(Er167)Di    77.77778
#> 126 10    CD38(Er167)Di     80.00000 10    CD38(Er167)Di    80.00000
#> 127 11    CD38(Er167)Di     81.81818 11    CD38(Er167)Di    81.81818
#> 128 12    CD38(Er167)Di     83.33333 12    CD38(Er167)Di    83.33333
#> 129  5    CD41(Lu175)Di     60.00000  5    CD41(Lu175)Di    60.00000
#> 130  6    CD41(Lu175)Di     66.66667  6    CD41(Lu175)Di    66.66667
#> 131  7    CD41(Lu175)Di     71.42857  7    CD41(Lu175)Di    71.42857
#> 132  8    CD41(Lu175)Di     75.00000  8    CD41(Lu175)Di    75.00000
#> 133  9    CD41(Lu175)Di     77.77778  9    CD41(Lu175)Di    77.77778
#> 134 10    CD41(Lu175)Di     80.00000 10    CD41(Lu175)Di    80.00000
#> 135 11    CD41(Lu175)Di     81.81818 11    CD41(Lu175)Di    81.81818
#> 136 12    CD41(Lu175)Di     83.33333 12    CD41(Lu175)Di    83.33333
#> 137  5    CD44(Er166)Di    100.00000  5    CD44(Er166)Di   100.00000
#> 138  6    CD44(Er166)Di    100.00000  6    CD44(Er166)Di   100.00000
#> 139  7    CD44(Er166)Di    100.00000  7    CD44(Er166)Di   100.00000
#> 140  8    CD44(Er166)Di    100.00000  8    CD44(Er166)Di   100.00000
#> 141  9    CD44(Er166)Di    100.00000  9    CD44(Er166)Di   100.00000
#> 142 10    CD44(Er166)Di    100.00000 10    CD44(Er166)Di   100.00000
#> 143 11    CD44(Er166)Di    100.00000 11    CD44(Er166)Di   100.00000
#> 144 12    CD44(Er166)Di    100.00000 12    CD44(Er166)Di   100.00000
#> 145  5    CD45(Sm154)Di    100.00000  5    CD45(Sm154)Di   100.00000
#> 146  6    CD45(Sm154)Di    100.00000  6    CD45(Sm154)Di   100.00000
#> 147  7    CD45(Sm154)Di    100.00000  7    CD45(Sm154)Di   100.00000
#> 148  8    CD45(Sm154)Di    100.00000  8    CD45(Sm154)Di   100.00000
#> 149  9    CD45(Sm154)Di    100.00000  9    CD45(Sm154)Di   100.00000
#> 150 10    CD45(Sm154)Di    100.00000 10    CD45(Sm154)Di   100.00000
#> 151 11    CD45(Sm154)Di    100.00000 11    CD45(Sm154)Di   100.00000
#> 152 12    CD45(Sm154)Di    100.00000 12    CD45(Sm154)Di   100.00000
#> 153  5  CD45RA(La139)Di     60.00000  5  CD45RA(La139)Di    60.00000
#> 154  6  CD45RA(La139)Di     66.66667  6  CD45RA(La139)Di    66.66667
#> 155  7  CD45RA(La139)Di     71.42857  7  CD45RA(La139)Di    71.42857
#> 156  8  CD45RA(La139)Di     75.00000  8  CD45RA(La139)Di    75.00000
#> 157  9  CD45RA(La139)Di     77.77778  9  CD45RA(La139)Di    77.77778
#> 158 10  CD45RA(La139)Di     80.00000 10  CD45RA(La139)Di    80.00000
#> 159 11  CD45RA(La139)Di     81.81818 11  CD45RA(La139)Di    81.81818
#> 160 12  CD45RA(La139)Di     83.33333 12  CD45RA(La139)Di    83.33333
#> 161  5    CD47(Gd160)Di    100.00000  5    CD47(Gd160)Di   100.00000
#> 162  6    CD47(Gd160)Di    100.00000  6    CD47(Gd160)Di   100.00000
#> 163  7    CD47(Gd160)Di    100.00000  7    CD47(Gd160)Di   100.00000
#> 164  8    CD47(Gd160)Di    100.00000  8    CD47(Gd160)Di   100.00000
#> 165  9    CD47(Gd160)Di    100.00000  9    CD47(Gd160)Di   100.00000
#> 166 10    CD47(Gd160)Di    100.00000 10    CD47(Gd160)Di   100.00000
#> 167 11    CD47(Gd160)Di    100.00000 11    CD47(Gd160)Di   100.00000
#> 168 12    CD47(Gd160)Di    100.00000 12    CD47(Gd160)Di   100.00000
#> 169  5   CD49d(Yb172)Di     80.00000  5   CD49d(Yb172)Di    80.00000
#> 170  6   CD49d(Yb172)Di     83.33333  6   CD49d(Yb172)Di    83.33333
#> 171  7   CD49d(Yb172)Di     85.71429  7   CD49d(Yb172)Di    85.71429
#> 172  8   CD49d(Yb172)Di     87.50000  8   CD49d(Yb172)Di    87.50000
#> 173  9   CD49d(Yb172)Di     88.88889  9   CD49d(Yb172)Di    88.88889
#> 174 10   CD49d(Yb172)Di     90.00000 10   CD49d(Yb172)Di    90.00000
#> 175 11   CD49d(Yb172)Di     90.90909 11   CD49d(Yb172)Di    90.90909
#> 176 12   CD49d(Yb172)Di     91.66667 12   CD49d(Yb172)Di    91.66667
#> 177  5    CD61(Tm169)Di     80.00000  5    CD61(Tm169)Di    80.00000
#> 178  6    CD61(Tm169)Di     83.33333  6    CD61(Tm169)Di    83.33333
#> 179  7    CD61(Tm169)Di     85.71429  7    CD61(Tm169)Di    85.71429
#> 180  8    CD61(Tm169)Di    100.00000  8    CD61(Tm169)Di   100.00000
#> 181  9    CD61(Tm169)Di    100.00000  9    CD61(Tm169)Di   100.00000
#> 182 10    CD61(Tm169)Di    100.00000 10    CD61(Tm169)Di   100.00000
#> 183 11    CD61(Tm169)Di    100.00000 11    CD61(Tm169)Di   100.00000
#> 184 12    CD61(Tm169)Di    100.00000 12    CD61(Tm169)Di   100.00000
#> 185  5    CD64(Yb176)Di     40.00000  5    CD64(Yb176)Di    40.00000
#> 186  6    CD64(Yb176)Di     50.00000  6    CD64(Yb176)Di    50.00000
#> 187  7    CD64(Yb176)Di     71.42857  7    CD64(Yb176)Di    71.42857
#> 188  8    CD64(Yb176)Di     62.50000  8    CD64(Yb176)Di    62.50000
#> 189  9    CD64(Yb176)Di     66.66667  9    CD64(Yb176)Di    66.66667
#> 190 10    CD64(Yb176)Di     70.00000 10    CD64(Yb176)Di    70.00000
#> 191 11    CD64(Yb176)Di     72.72727 11    CD64(Yb176)Di    72.72727
#> 192 12    CD64(Yb176)Di     66.66667 12    CD64(Yb176)Di    66.66667
#> 193  5   CD117(Yb171)Di    100.00000  5   CD117(Yb171)Di   100.00000
#> 194  6   CD117(Yb171)Di    100.00000  6   CD117(Yb171)Di   100.00000
#> 195  7   CD117(Yb171)Di    100.00000  7   CD117(Yb171)Di   100.00000
#> 196  8   CD117(Yb171)Di    100.00000  8   CD117(Yb171)Di   100.00000
#> 197  9   CD117(Yb171)Di    100.00000  9   CD117(Yb171)Di   100.00000
#> 198 10   CD117(Yb171)Di    100.00000 10   CD117(Yb171)Di   100.00000
#> 199 11   CD117(Yb171)Di    100.00000 11   CD117(Yb171)Di   100.00000
#> 200 12   CD117(Yb171)Di    100.00000 12   CD117(Yb171)Di   100.00000
#> 201  5   CD123(Eu151)Di      0.00000  5   CD123(Eu151)Di     0.00000
#> 202  6   CD123(Eu151)Di      0.00000  6   CD123(Eu151)Di     0.00000
#> 203  7   CD123(Eu151)Di     14.28571  7   CD123(Eu151)Di    14.28571
#> 204  8   CD123(Eu151)Di     12.50000  8   CD123(Eu151)Di    12.50000
#> 205  9   CD123(Eu151)Di     22.22222  9   CD123(Eu151)Di    22.22222
#> 206 10   CD123(Eu151)Di     40.00000 10   CD123(Eu151)Di    40.00000
#> 207 11   CD123(Eu151)Di     54.54545 11   CD123(Eu151)Di    54.54545
#> 208 12   CD123(Eu151)Di     50.00000 12   CD123(Eu151)Di    50.00000
#> 209  5   CD133(Pr141)Di     20.00000  5   CD133(Pr141)Di    20.00000
#> 210  6   CD133(Pr141)Di     33.33333  6   CD133(Pr141)Di    33.33333
#> 211  7   CD133(Pr141)Di     57.14286  7   CD133(Pr141)Di    57.14286
#> 212  8   CD133(Pr141)Di     50.00000  8   CD133(Pr141)Di    50.00000
#> 213  9   CD133(Pr141)Di     55.55556  9   CD133(Pr141)Di    55.55556
#> 214 10   CD133(Pr141)Di     50.00000 10   CD133(Pr141)Di    50.00000
#> 215 11   CD133(Pr141)Di     54.54545 11   CD133(Pr141)Di    54.54545
#> 216 12   CD133(Pr141)Di     58.33333 12   CD133(Pr141)Di    58.33333
#> 217  5 CD235ab(Sm152)Di     60.00000  5 CD235ab(Sm152)Di    60.00000
#> 218  6 CD235ab(Sm152)Di     66.66667  6 CD235ab(Sm152)Di    66.66667
#> 219  7 CD235ab(Sm152)Di     71.42857  7 CD235ab(Sm152)Di    71.42857
#> 220  8 CD235ab(Sm152)Di     75.00000  8 CD235ab(Sm152)Di    75.00000
#> 221  9 CD235ab(Sm152)Di     77.77778  9 CD235ab(Sm152)Di    77.77778
#> 222 10 CD235ab(Sm152)Di     80.00000 10 CD235ab(Sm152)Di    80.00000
#> 223 11 CD235ab(Sm152)Di     81.81818 11 CD235ab(Sm152)Di    81.81818
#> 224 12 CD235ab(Sm152)Di     83.33333 12 CD235ab(Sm152)Di    83.33333
#> 225  5   CD321(Eu153)Di    100.00000  5   CD321(Eu153)Di   100.00000
#> 226  6   CD321(Eu153)Di    100.00000  6   CD321(Eu153)Di   100.00000
#> 227  7   CD321(Eu153)Di    100.00000  7   CD321(Eu153)Di   100.00000
#> 228  8   CD321(Eu153)Di    100.00000  8   CD321(Eu153)Di   100.00000
#> 229  9   CD321(Eu153)Di    100.00000  9   CD321(Eu153)Di   100.00000
#> 230 10   CD321(Eu153)Di    100.00000 10   CD321(Eu153)Di   100.00000
#> 231 11   CD321(Eu153)Di    100.00000 11   CD321(Eu153)Di   100.00000
#> 232 12   CD321(Eu153)Di    100.00000 12   CD321(Eu153)Di   100.00000
#> 233  5   CXCR4(Sm149)Di     60.00000  5   CXCR4(Sm149)Di    60.00000
#> 234  6   CXCR4(Sm149)Di     50.00000  6   CXCR4(Sm149)Di    50.00000
#> 235  7   CXCR4(Sm149)Di     57.14286  7   CXCR4(Sm149)Di    57.14286
#> 236  8   CXCR4(Sm149)Di     62.50000  8   CXCR4(Sm149)Di    62.50000
#> 237  9   CXCR4(Sm149)Di     66.66667  9   CXCR4(Sm149)Di    66.66667
#> 238 10   CXCR4(Sm149)Di     70.00000 10   CXCR4(Sm149)Di    70.00000
#> 239 11   CXCR4(Sm149)Di     72.72727 11   CXCR4(Sm149)Di    72.72727
#> 240 12   CXCR4(Sm149)Di     83.33333 12   CXCR4(Sm149)Di    83.33333
#> 241  5    Flt3(Nd150)Di     20.00000  5    Flt3(Nd150)Di    20.00000
#> 242  6    Flt3(Nd150)Di     33.33333  6    Flt3(Nd150)Di    33.33333
#> 243  7    Flt3(Nd150)Di     28.57143  7    Flt3(Nd150)Di    28.57143
#> 244  8    Flt3(Nd150)Di     25.00000  8    Flt3(Nd150)Di    25.00000
#> 245  9    Flt3(Nd150)Di     22.22222  9    Flt3(Nd150)Di    22.22222
#> 246 10    Flt3(Nd150)Di     30.00000 10    Flt3(Nd150)Di    30.00000
#> 247 11    Flt3(Nd150)Di     45.45455 11    Flt3(Nd150)Di    45.45455
#> 248 12    Flt3(Nd150)Di     50.00000 12    Flt3(Nd150)Di    50.00000
#> 249  5  HLA-DR(Yb174)Di     40.00000  5  HLA-DR(Yb174)Di    40.00000
#> 250  6  HLA-DR(Yb174)Di     50.00000  6  HLA-DR(Yb174)Di    50.00000
#> 251  7  HLA-DR(Yb174)Di     42.85714  7  HLA-DR(Yb174)Di    42.85714
#> 252  8  HLA-DR(Yb174)Di     37.50000  8  HLA-DR(Yb174)Di    37.50000
#> 253  9  HLA-DR(Yb174)Di     33.33333  9  HLA-DR(Yb174)Di    33.33333
#> 254 10  HLA-DR(Yb174)Di     40.00000 10  HLA-DR(Yb174)Di    40.00000
#> 255 11  HLA-DR(Yb174)Di     45.45455 11  HLA-DR(Yb174)Di    45.45455
#> 256 12  HLA-DR(Yb174)Di     41.66667 12  HLA-DR(Yb174)Di    41.66667
#> 
#> $plot

#> 

```
