/* Re-emit a numeric TSV with every field as a C99 hex float (issue #13).
 *
 * R's decimal string-to-double conversion is not correctly rounded. `as.numeric`,
 * `scan` and `read.delim` all share it, and on a random sample it returns a
 * double one ULP away from the correctly rounded value for about 18% of inputs:
 *
 *     "26.550866314209998"   R gives 0x1.a8d0593240001p+4
 *                            correct  0x1.a8d059324p+4
 *
 * So a golden file written by C++ at full `%.17g` precision does NOT survive
 * being read into R -- roughly a fifth of its values arrive one ULP off. Any
 * comparison that reads C++ output through R while computing the other side
 * inside R is therefore measuring R's parser, not the models.
 *
 * This filter does the parse with the C library's strtod, which IS correctly
 * rounded, and prints the result as hex. R reads hex exactly (verified: 4000 of
 * 4000 values round-trip, against 3265 of 4000 for %.17g), so the values reach R
 * unharmed.
 *
 *   cc -O2 -o tsv_to_hex tests/validate/tsv_to_hex.c
 *   tsv_to_hex < golden/operating_points.tsv > golden_hex.tsv
 *
 * Fields that do not parse as a number are passed through unchanged, so a header
 * row needs no special handling.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE 65536

static void emit_field(const char *s) {
  char *end = NULL;
  double v = strtod(s, &end);
  /* A number only if strtod consumed the whole field and the field was not
     empty. Anything else -- a column name, a stray token -- goes out as-is. */
  if (end != s && *end == '\0') {
    printf("%a", v);
  } else {
    fputs(s, stdout);
  }
}

int main(void) {
  char line[MAX_LINE];
  while (fgets(line, sizeof line, stdin)) {
    size_t n = strlen(line);
    int had_newline = (n > 0 && line[n - 1] == '\n');
    if (had_newline) line[--n] = '\0';
    if (n > 0 && line[n - 1] == '\r') line[--n] = '\0';

    const char *start = line;
    for (;;) {
      char *tab = strchr(start, '\t');
      if (tab) *tab = '\0';
      emit_field(start);
      if (!tab) break;
      putchar('\t');
      start = tab + 1;
    }
    if (had_newline) putchar('\n');
  }
  return 0;
}
