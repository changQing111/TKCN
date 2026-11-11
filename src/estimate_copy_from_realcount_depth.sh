#!/usr/bin/env bash
set -euo pipefail

target_inputs=()
depth_file=""
out_file=""
target_mad_z=3
target_method="mad_filtered_median"
copy_half_window=""
copy_call_method=""
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python_bin=${PYTHON:-/data/changq/miniforge3/bin/python}
nb_helper="$script_dir/estimate_depth_from_realcount_nb.py"

usage() {
  cat <<EOF
Usage:
  bash $0 --target FILE_OR_DIR --depth SINGLE_COPY_DEPTH.tsv [options]

Estimate target-gene copy number from target-specific kmer realcount files
and sequencing depth estimated from a single-copy marker.

Formula:
  default:
    copy_estimate = median_mad_filtered(target_specific_kmer_counts) / single_copy_depth
  with --target-method nb:
    copy_estimate = target_nb_mean_count / single_copy_depth

Required:
  -t, --target PATH   Target realcount TSV file or directory. Can be repeated.
                      Directory inputs use *.realcount.tsv and *_inter.k31.tsv.
  -d, --depth FILE    Depth TSV with depth_estimate or nb_depth_estimate.

Options:
  -o, --out FILE      Write result TSV to FILE. Default: stdout.
  --target-method METHOD
                      Target kmer summary method: mad_filtered_median or nb.
                      Default: $target_method.
  --target-mad-z NUM  MAD filter width for target kmer counts. Keeps counts
                      within median +/- NUM * 1.4826 * MAD. Default: $target_mad_z.
  --copy-half-window NUM
                      Add copy-call columns and mark estimates within NUM of
                      a half-integer as ambiguous. Example: 0.1 marks
                      2.40-2.60 as candidate copies 2,3.
  --copy-call-method METHOD
                      Copy calling method for extra copy_call columns:
                      half_window or nb_likelihood. nb_likelihood requires
                      --target-method nb and forces an integer by comparing
                      adjacent integer-copy NB/Poisson likelihoods.
  -h, --help          Show this help.

Examples:
  bash $0 --target TaFT-B1_num/85_inter.k31.tsv --depth single_copy_depth_estimate.tsv
  bash $0 --target multi_copy_gene_kmer_num/TaFT-B1_copy_num \\
    --depth single_copy_gene_kmer_num/single_copy_depth_estimate.tsv \\
    --out TaFT-B1_copy_estimate.tsv
  bash $0 --target multi_copy_gene_kmer_num/TaFT-B1_copy_num \\
    --depth single_copy_gene_kmer_num/single_copy_depth_nb_joint.tsv \\
    --target-method nb \\
    --out TaFT-B1_copy_estimate_target_nb.tsv
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -t|--target)
      [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a file or directory" >&2; exit 1; }
      target_inputs+=("$2")
      shift 2
      ;;
    --target=*)
      target_inputs+=("${1#*=}")
      shift
      ;;
    -d|--depth)
      [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a depth TSV" >&2; exit 1; }
      depth_file=$2
      shift 2
      ;;
    --depth=*)
      depth_file=${1#*=}
      shift
      ;;
    -o|--out)
      [[ $# -ge 2 ]] || { echo "ERROR: $1 requires an output file" >&2; exit 1; }
      out_file=$2
      shift 2
      ;;
    --out=*)
      out_file=${1#*=}
      shift
      ;;
    --target-mad-z|--target_mad_z)
      [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a non-negative number" >&2; exit 1; }
      target_mad_z=$2
      shift 2
      ;;
    --target-mad-z=*|--target_mad_z=*)
      target_mad_z=${1#*=}
      shift
      ;;
    --target-method|--target_method)
      [[ $# -ge 2 ]] || { echo "ERROR: $1 requires mad_filtered_median or nb" >&2; exit 1; }
      target_method=$2
      shift 2
      ;;
    --target-method=*|--target_method=*)
      target_method=${1#*=}
      shift
      ;;
    --copy-half-window|--copy_half_window)
      [[ $# -ge 2 ]] || { echo "ERROR: $1 requires a number from 0 to 0.5" >&2; exit 1; }
      copy_half_window=$2
      shift 2
      ;;
    --copy-half-window=*|--copy_half_window=*)
      copy_half_window=${1#*=}
      shift
      ;;
    --copy-call-method|--copy_call_method)
      [[ $# -ge 2 ]] || { echo "ERROR: $1 requires half_window or nb_likelihood" >&2; exit 1; }
      copy_call_method=$2
      shift 2
      ;;
    --copy-call-method=*|--copy_call_method=*)
      copy_call_method=${1#*=}
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      echo "ERROR: unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
    *)
      echo "ERROR: unexpected positional argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ ${#target_inputs[@]} -eq 0 ]]; then
  echo "ERROR: --target is required" >&2
  usage >&2
  exit 1
fi

if [[ -z "$depth_file" ]]; then
  echo "ERROR: --depth is required" >&2
  usage >&2
  exit 1
fi

if [[ ! -f "$depth_file" ]]; then
  echo "ERROR: depth file not found: $depth_file" >&2
  exit 1
fi

if ! awk -v x="$target_mad_z" 'BEGIN { exit !(x ~ /^[0-9]+([.][0-9]+)?$/ && x >= 0) }'; then
  echo "ERROR: --target-mad-z must be a non-negative number: $target_mad_z" >&2
  exit 1
fi

if [[ -n "$copy_half_window" ]] && ! awk -v x="$copy_half_window" 'BEGIN { exit !(x ~ /^[0-9]+([.][0-9]+)?$/ && x >= 0 && x <= 0.5) }'; then
  echo "ERROR: --copy-half-window must be a number from 0 to 0.5: $copy_half_window" >&2
  exit 1
fi

if [[ -n "$copy_half_window" && -z "$copy_call_method" ]]; then
  copy_call_method="half_window"
fi

case "$target_method" in
  mad_filtered_median|nb)
    ;;
  *)
    echo "ERROR: --target-method must be mad_filtered_median or nb: $target_method" >&2
    exit 1
    ;;
esac

case "$copy_call_method" in
  ""|half_window|nb_likelihood)
    ;;
  *)
    echo "ERROR: --copy-call-method must be half_window or nb_likelihood: $copy_call_method" >&2
    exit 1
    ;;
esac

if [[ "$copy_call_method" == "half_window" && -z "$copy_half_window" ]]; then
  echo "ERROR: --copy-call-method half_window requires --copy-half-window" >&2
  exit 1
fi

if [[ "$copy_call_method" == "nb_likelihood" && "$target_method" != "nb" ]]; then
  echo "ERROR: --copy-call-method nb_likelihood requires --target-method nb" >&2
  exit 1
fi

expanded_targets=()
for input in "${target_inputs[@]}"; do
  if [[ -d "$input" ]]; then
    while IFS= read -r path; do
      expanded_targets+=("$path")
    done < <(
      find -L "$input" -maxdepth 1 -type f \
        \( -name "*.realcount.tsv" -o -name "*_inter.k31.tsv" \)
    )
  elif [[ -f "$input" ]]; then
    expanded_targets+=("$input")
  else
    echo "ERROR: target input not found: $input" >&2
    exit 1
  fi
done

target_files=()
if [[ ${#expanded_targets[@]} -gt 0 ]]; then
  mapfile -t target_files < <(printf '%s\n' "${expanded_targets[@]}" | sort -u)
fi

if [[ ${#target_files[@]} -eq 0 ]]; then
  echo "ERROR: no *.realcount.tsv or *_inter.k31.tsv target files found" >&2
  exit 1
fi

sample_from_path() {
  local path=$1
  local name sample

  name=$(basename "$path")
  sample=$name
  sample=${sample%.tsv}
  sample=${sample%_inter.k31}
  sample=${sample%.realcount}
  sample=${sample%.panel}

  printf '%s' "$sample"
}

target_summary_tmp=$(mktemp)
target_nb_tmp=$(mktemp)
result_tmp=$(mktemp)
cleanup() {
  rm -f "$target_summary_tmp" "$target_nb_tmp" "$result_tmp"
}
trap cleanup EXIT

if [[ "$target_method" == "mad_filtered_median" ]]; then
  for target_file in "${target_files[@]}"; do
    sample=$(sample_from_path "$target_file")
    awk \
      -v sample="$sample" \
      -v file="$target_file" \
      -v target_mad_z="$target_mad_z" \
      '
      function format_number(x) {
        return sprintf("%.10g", x)
      }

      function median_range(arr, start, stop, middle) {
        middle = int((start + stop) / 2)
        if (((stop - start + 1) % 2) == 1) {
          return arr[middle]
        }
        return (arr[middle] + arr[middle + 1]) / 2
      }

      function fail(message) {
        printf "ERROR: %s in %s at line %d\n", message, file, NR > "/dev/stderr"
        exit 1
      }

      /^[[:space:]]*$/ {
        next
      }

      NF < 2 {
        fail("expected at least 2 columns")
      }

      $2 !~ /^[0-9]+([.][0-9]+)?$/ {
        fail("non-numeric count")
      }

      {
        n++
        counts[n] = $2 + 0
      }

      END {
        if (n == 0) {
          printf "ERROR: no kmer counts in %s\n", file > "/dev/stderr"
          exit 1
        }

        asort(counts)
        median_raw = median_range(counts, 1, n)

        for (i = 1; i <= n; i++) {
          diff[i] = counts[i] - median_raw
          if (diff[i] < 0) {
            diff[i] = -diff[i]
          }
        }
        asort(diff)
        mad = median_range(diff, 1, n)
        mad_radius = target_mad_z * 1.4826 * mad
        mad_lower = median_raw - mad_radius
        mad_upper = median_raw + mad_radius

        filtered_n = 0
        for (i = 1; i <= n; i++) {
          if (counts[i] >= mad_lower && counts[i] <= mad_upper) {
            filtered_n++
            filtered[filtered_n] = counts[i]
          }
        }

        if (filtered_n == 0) {
          filtered_n = n
          for (i = 1; i <= n; i++) {
            filtered[i] = counts[i]
          }
        }

        median_mad_filtered = median_range(filtered, 1, filtered_n)

        printf "%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\n", \
          sample, file, n, format_number(median_mad_filtered), filtered_n, \
          format_number(mad), format_number(mad_lower), format_number(mad_upper)
      }
      ' "$target_file" >> "$target_summary_tmp"
  done
else
  "$python_bin" "$nb_helper" marker --copy-number 1 "${target_files[@]}" > "$target_nb_tmp"
fi

write_output() {
  if [[ "$target_method" == "mad_filtered_median" ]]; then
    if [[ "$copy_call_method" == "half_window" ]]; then
      printf 'sample\ttarget_file\tdepth_file\ttarget_kmer_count\ttarget_median_mad_filtered\ttarget_mad_filtered_kmer_count\ttarget_mad\ttarget_mad_lower\ttarget_mad_upper\tsingle_copy_depth\tcopy_estimate\tcopy_estimate_round\tcopy_call\tcopy_call_status\tcopy_call_candidates\tcopy_fractional_part\tcopy_distance_to_half_integer\n' > "$result_tmp"
    else
      printf 'sample\ttarget_file\tdepth_file\ttarget_kmer_count\ttarget_median_mad_filtered\ttarget_mad_filtered_kmer_count\ttarget_mad\ttarget_mad_lower\ttarget_mad_upper\tsingle_copy_depth\tcopy_estimate\tcopy_estimate_round\n' > "$result_tmp"
    fi

    awk \
      -v depth_file="$depth_file" \
      -v result_file="$result_tmp" \
      -v copy_half_window="$copy_half_window" \
      -v copy_call_method="$copy_call_method" \
      '
        function format_number(x) {
          return sprintf("%.10g", x)
        }

        function round_nearest(x) {
          return int(x + 0.5)
        }

        function append_copy_call(copy_estimate, rounded, lower, fractional, distance, candidate_a, candidate_b) {
          lower = int(copy_estimate)
          fractional = copy_estimate - lower
          distance = fractional - 0.5
          if (distance < 0) {
            distance = -distance
          }
          if (distance <= copy_half_window + 0) {
            candidate_a = lower
            candidate_b = lower + 1
            return OFS "NA" OFS "ambiguous_half_integer" OFS candidate_a "," candidate_b OFS format_number(fractional) OFS format_number(distance)
          }
          return OFS rounded OFS "ok" OFS rounded OFS format_number(fractional) OFS format_number(distance)
        }

        BEGIN {
          FS = OFS = "\t"
        }

        NR == FNR {
          if (FNR == 1) {
            for (i = 1; i <= NF; i++) {
              header[$i] = i
            }
            if (!("sample" in header)) {
              printf "ERROR: depth file must contain a sample column: %s\n", depth_file > "/dev/stderr"
              exit 1
            }
            if ("depth_estimate" in header) {
              depth_column = "depth_estimate"
            } else if ("nb_depth_estimate" in header) {
              depth_column = "nb_depth_estimate"
            } else {
              printf "ERROR: depth file must contain depth_estimate or nb_depth_estimate column: %s\n", depth_file > "/dev/stderr"
              exit 1
            }
            next
          }

          sample = $header["sample"]
          depth = $header[depth_column]
          if (depth !~ /^[0-9]+([.][0-9]+)?$/ || depth <= 0) {
            printf "ERROR: invalid %s for sample %s in %s\n", depth_column, sample, depth_file > "/dev/stderr"
            exit 1
          }
          depths[sample] = depth + 0
          next
        }

        {
          sample = $1
          target_file = $2
          target_kmer_count = $3
          target_median_mad_filtered = $4 + 0
          target_mad_filtered_kmer_count = $5
          target_mad = $6 + 0
          target_mad_lower = $7 + 0
          target_mad_upper = $8 + 0

          if (!(sample in depths)) {
            printf "WARNING: skipped sample without depth: %s\n", sample > "/dev/stderr"
            next
          }

          copy_estimate = target_median_mad_filtered / depths[sample]
          rounded = round_nearest(copy_estimate)
          line = sample OFS target_file OFS depth_file OFS target_kmer_count OFS \
            format_number(target_median_mad_filtered) OFS target_mad_filtered_kmer_count OFS \
            format_number(target_mad) OFS format_number(target_mad_lower) OFS format_number(target_mad_upper) OFS \
            format_number(depths[sample]) OFS \
            format_number(copy_estimate) OFS rounded
          if (copy_call_method == "half_window") {
            line = line append_copy_call(copy_estimate, rounded)
          }
          print line >> result_file
        }
      ' "$depth_file" "$target_summary_tmp"
  else
    if [[ "$copy_call_method" == "half_window" ]]; then
      printf 'sample\ttarget_file\tdepth_file\ttarget_kmer_count\ttarget_nb_mean_count\ttarget_nb_var_count\ttarget_nb_size\ttarget_nb_loglik\ttarget_nb_status\tsingle_copy_depth\tcopy_estimate\tcopy_estimate_round\tcopy_call\tcopy_call_status\tcopy_call_candidates\tcopy_fractional_part\tcopy_distance_to_half_integer\n' > "$result_tmp"
    elif [[ "$copy_call_method" == "nb_likelihood" ]]; then
      printf 'sample\ttarget_file\tdepth_file\ttarget_kmer_count\ttarget_nb_mean_count\ttarget_nb_var_count\ttarget_nb_size\ttarget_nb_loglik\ttarget_nb_status\tsingle_copy_depth\tcopy_estimate\tcopy_estimate_round\tcopy_call\tcopy_call_status\tcopy_call_candidates\tcopy_fractional_part\tcopy_distance_to_half_integer\tcopy_call_logliks\tcopy_call_loglik_delta\n' > "$result_tmp"
    else
      printf 'sample\ttarget_file\tdepth_file\ttarget_kmer_count\ttarget_nb_mean_count\ttarget_nb_var_count\ttarget_nb_size\ttarget_nb_loglik\ttarget_nb_status\tsingle_copy_depth\tcopy_estimate\tcopy_estimate_round\n' > "$result_tmp"
    fi

    awk \
      -v depth_file="$depth_file" \
      -v result_file="$result_tmp" \
      -v copy_half_window="$copy_half_window" \
      -v copy_call_method="$copy_call_method" \
      '
        function format_number(x) {
          return sprintf("%.10g", x)
        }

        function round_nearest(x) {
          return int(x + 0.5)
        }

        function append_copy_call(copy_estimate, rounded, lower, fractional, distance, candidate_a, candidate_b) {
          lower = int(copy_estimate)
          fractional = copy_estimate - lower
          distance = fractional - 0.5
          if (distance < 0) {
            distance = -distance
          }
          if (distance <= copy_half_window + 0) {
            candidate_a = lower
            candidate_b = lower + 1
            return OFS "NA" OFS "ambiguous_half_integer" OFS candidate_a "," candidate_b OFS format_number(fractional) OFS format_number(distance)
          }
          return OFS rounded OFS "ok" OFS rounded OFS format_number(fractional) OFS format_number(distance)
        }

        function is_numeric(x) {
          return x ~ /^[-+]?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$/ || x ~ /^[-+]?[.][0-9]+([eE][-+]?[0-9]+)?$/
        }

        function copy_model_loglik(candidate, depth, kmer_count, sum_count, nb_size, nb_status, nb_loglik_constant, mu, size, constant) {
          mu = candidate * depth
          if (candidate == 0) {
            if (sum_count == 0) {
              return 0
            }
            return -1e300
          }
          if (mu <= 0 || !is_numeric(nb_loglik_constant)) {
            return ""
          }
          constant = nb_loglik_constant + 0
          if (nb_status == "ok" && is_numeric(nb_size)) {
            size = nb_size + 0
            if (size <= 0) {
              return ""
            }
            return constant + sum_count * log(mu / (size + mu)) + kmer_count * size * log(size / (size + mu))
          }
          if (nb_status == "poisson_limit") {
            return constant + sum_count * log(mu) - kmer_count * mu
          }
          return ""
        }

        function append_nb_likelihood_copy_call(copy_estimate, rounded, depth, kmer_count, sum_count, nb_size, nb_status, nb_loglik_constant, lower, upper, fractional, distance, lower_ll, upper_ll, best, delta) {
          lower = int(copy_estimate)
          fractional = copy_estimate - lower
          upper = lower + 1
          distance = fractional - 0.5
          if (distance < 0) {
            distance = -distance
          }
          if (fractional == 0) {
            upper = lower
          }

          lower_ll = copy_model_loglik(lower, depth, kmer_count, sum_count, nb_size, nb_status, nb_loglik_constant)
          upper_ll = copy_model_loglik(upper, depth, kmer_count, sum_count, nb_size, nb_status, nb_loglik_constant)
          if (lower_ll == "" || upper_ll == "") {
            return OFS rounded OFS "nb_likelihood_unavailable" OFS lower "," upper OFS format_number(fractional) OFS format_number(distance) OFS "NA" OFS "NA"
          }
          if ((lower_ll + 0) >= (upper_ll + 0)) {
            best = lower
            delta = (lower_ll + 0) - (upper_ll + 0)
          } else {
            best = upper
            delta = (upper_ll + 0) - (lower_ll + 0)
          }
          return OFS best OFS "nb_likelihood" OFS lower "," upper OFS format_number(fractional) OFS format_number(distance) OFS format_number(lower_ll) "," format_number(upper_ll) OFS format_number(delta)
        }

        BEGIN {
          FS = OFS = "\t"
        }

        NR == FNR {
          if (FNR == 1) {
            for (i = 1; i <= NF; i++) {
              depth_header[$i] = i
            }
            if (!("sample" in depth_header)) {
              printf "ERROR: depth file must contain a sample column: %s\n", depth_file > "/dev/stderr"
              exit 1
            }
            if ("depth_estimate" in depth_header) {
              depth_column = "depth_estimate"
            } else if ("nb_depth_estimate" in depth_header) {
              depth_column = "nb_depth_estimate"
            } else {
              printf "ERROR: depth file must contain depth_estimate or nb_depth_estimate column: %s\n", depth_file > "/dev/stderr"
              exit 1
            }
            next
          }

          sample = $depth_header["sample"]
          depth = $depth_header[depth_column]
          if (depth !~ /^[0-9]+([.][0-9]+)?$/ || depth <= 0) {
            printf "ERROR: invalid %s for sample %s in %s\n", depth_column, sample, depth_file > "/dev/stderr"
            exit 1
          }
          depths[sample] = depth + 0
          next
        }

        FNR == 1 {
          for (i = 1; i <= NF; i++) {
            target_header[$i] = i
          }
          required[1] = "sample"
          required[2] = "input_file"
          required[3] = "kmer_count"
          required[4] = "mean_count"
          required[5] = "var_count"
          required[6] = "nb_size"
          required[7] = "nb_loglik"
          required[8] = "nb_status"
          required[9] = "sum_count"
          required[10] = "nb_loglik_constant"
          for (i = 1; i <= 10; i++) {
            if (!(required[i] in target_header)) {
              printf "ERROR: target NB table missing %s column\n", required[i] > "/dev/stderr"
              exit 1
            }
          }
          next
        }

        {
          sample = $target_header["sample"]
          target_file = $target_header["input_file"]
          target_kmer_count = $target_header["kmer_count"]
          target_sum_count = $target_header["sum_count"]
          target_nb_mean_count = $target_header["mean_count"]
          target_nb_var_count = $target_header["var_count"]
          target_nb_size = $target_header["nb_size"]
          target_nb_loglik = $target_header["nb_loglik"]
          target_nb_status = $target_header["nb_status"]
          target_nb_loglik_constant = $target_header["nb_loglik_constant"]

          if (!(sample in depths)) {
            printf "WARNING: skipped sample without depth: %s\n", sample > "/dev/stderr"
            next
          }
          if (target_nb_mean_count !~ /^[0-9]+([.][0-9]+)?$/) {
            printf "ERROR: invalid target NB mean_count for sample %s\n", sample > "/dev/stderr"
            exit 1
          }
          if (target_sum_count !~ /^[0-9]+([.][0-9]+)?$/ || target_kmer_count !~ /^[0-9]+$/) {
            printf "ERROR: invalid target NB count fields for sample %s\n", sample > "/dev/stderr"
            exit 1
          }

          copy_estimate = (target_nb_mean_count + 0) / depths[sample]
          rounded = round_nearest(copy_estimate)
          line = sample OFS target_file OFS depth_file OFS target_kmer_count OFS \
            target_nb_mean_count OFS target_nb_var_count OFS target_nb_size OFS target_nb_loglik OFS target_nb_status OFS \
            format_number(depths[sample]) OFS \
            format_number(copy_estimate) OFS rounded
          if (copy_call_method == "half_window") {
            line = line append_copy_call(copy_estimate, rounded)
          } else if (copy_call_method == "nb_likelihood") {
            line = line append_nb_likelihood_copy_call(copy_estimate, rounded, depths[sample], target_kmer_count + 0, target_sum_count + 0, target_nb_size, target_nb_status, target_nb_loglik_constant)
          }
          print line >> result_file
        }
      ' "$depth_file" "$target_nb_tmp"
  fi

  cat "$result_tmp"
}

if [[ -n "$out_file" ]]; then
  mkdir -p "$(dirname "$out_file")"
  write_output > "$out_file"
  written_count=$(awk 'NR > 1 {count++} END {print count + 0}' "$out_file")
  echo "Wrote $written_count copy estimates to $out_file" >&2
else
  write_output
fi
