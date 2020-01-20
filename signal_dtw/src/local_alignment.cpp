#include <local_alignment.h>
#include <matrix_band.h>
#include <limits>
#include <queue>
#include <fstream>
#include <cmath>
using namespace std;

double penalty(double value1, double value2) {
  double diff = value1 - value2;
  return diff * diff;
}

double prob_sum(double a, double b) {
  if (a < b) swap(a, b);
  return log(1 + exp(b - a)) + a;
}

double prob_rest(double a) {
  return log(1.0 - exp(a));
}

constexpr double kPSame = 0.99;
constexpr double kPDifferent = 1 - kPSame;
constexpr double kPSameMinusDifferent = 2*kPSame - 1;
constexpr double kLogPSame = log(kPSame);
constexpr double kLogPDifferent = log(kPDifferent);
constexpr double kLogPSameMinusDifferent = log(2*kPSame - 1);

double penalty(const vector<float> &predictions1, const vector<float> &predictions2) {
  double product = 0;
  for (unsigned i = 0; i < predictions1.size(); i++) {
    product += predictions1[i] * predictions2[i];
  }
  return -log(product * kPSameMinusDifferent + kPDifferent);
}

const int kMaxEventLength = 15;

vector<int> split_to_events(const vector<double>& signal, double event_threshold) {
  vector<int> result(1, 0);
  
  
  double max_val = numeric_limits<double>::lowest(), min_val = numeric_limits<double>::max();
  for (int i = 0; i < static_cast<int>(signal.size()); i++) {
    max_val = max(max_val, signal[i]);
    min_val = min(min_val, signal[i]);
    if (max_val - min_val > event_threshold || i >= result.back() + kMaxEventLength) {
      result.push_back(i);
      min_val = signal[i];
      max_val = signal[i];
    }
  }
  result.push_back(signal.size());
  return result;
}

vector<vector<pair<int, int>>> local_alignment(vector<double> signal,
                                               vector<vector<float>> chiron_predictions,
                                               int min_events_distance,
                                               int max_events_distance,
                                               double score_for_moving,
                                               int max_speed_ratio,
                                               double event_threshold,
                                               int min_lookahead,
                                               string log_filename) {
  fprintf(stderr, "starting local alignment by signal segmentation\n");
  vector<int> event_starts = split_to_events(signal, event_threshold * 2.0 /*this is a Bulgarian constant*/);
  fprintf(stderr, "signal segmentated, computing row starts and ends\n");
  int events_count = static_cast<int>(event_starts.size()) - 1;
  
  vector<int> row_starts(signal.size() + 1), row_ends(signal.size() + 1);
  vector<int> sum_appears(signal.size() * 2 + 1, -1);
  vector<int> col_appears(signal.size() + 1, -1);
  int next_sum = 0;
  int next_col = 0;
  int current_event_index = 0;
  for (int row = 0; row <= static_cast<int>(signal.size()); row++) {
    row_starts[row] = event_starts[min(current_event_index + min_events_distance, events_count)];
    row_ends[row] = event_starts[min(current_event_index + max_events_distance, events_count)] + 1;
    for (int sum = max(row + row_starts[row], next_sum); sum < row + row_ends[row]; sum++) {
      sum_appears[sum] = row;
      next_sum = sum+1;
    }
    for (int column = max(row_starts[row], next_col); column < row_ends[row]; column++) {
      col_appears[column] = row;
      next_col = column+1;
    }
    if (event_starts[current_event_index+1] == row) current_event_index++;
  }
  
  fprintf(stderr, "rows computed, computing steps back\n");
  vector<int> max_steps_back(signal.size()+1, 0);
  for (int i = 1; i <= static_cast<int>(signal.size()); i++) {
    double min_val = signal[i-1], max_val = signal[i-1];
    int steps = 1;
    for (int j = i-2; j >= row_starts[i] || (col_appears[i] != -1 && j >= col_appears[i]); j--) {
      min_val = min(min_val, signal[j]);
      max_val = max(max_val, signal[j]);
      if (max_val - min_val <= event_threshold) {
        steps++;
      }
      else break;
    }
    max_steps_back[i] = max(min(max_speed_ratio, i), steps);
  }
  
  MatrixBand<double> score(row_starts, row_ends, 0, numeric_limits<double>::lowest());
  MatrixBand<pair<int, int>> come_from(row_starts, row_ends, {-1, -1}, {-1, -1});
  
  vector<float> filler_prediction = {0, 0, 0, 0, 1};
  
  fprintf(stderr, "steps back computed, running dtw\n");
  int min_index_sum = row_starts[0];
  int max_index_sum = row_ends.back()-1 + signal.size();
  for (int window_start = min_index_sum; window_start + min_lookahead <= max_index_sum || window_start == min_index_sum; window_start += min_lookahead) {
    fprintf(stderr, "window start: %d / %d\r", window_start, max_index_sum);
    for (int index_sum = window_start; index_sum <= window_start + 2*min_lookahead; index_sum++) {
      bool update_come_from = (index_sum >= window_start + min_lookahead) || (window_start == min_index_sum);
      
      if (index_sum >= static_cast<int>(sum_appears.size()) || sum_appears[index_sum] == -1) continue;
      for (int row = sum_appears[index_sum]; row <= static_cast<int>(signal.size()) && row + row_starts[row] <= index_sum; row++) {
        int column = index_sum - row;
        
        double &current_cell = score[row][column];
        current_cell = 0;
        
        if (column - 1 >= row_starts[row]) {
          double step_score = score_for_moving; // for moving by one column
        
          int from_col = column - 1;
          for (int row_step = 1; row_step <= max_steps_back[row] && row_ends[row - row_step] > column; row_step++) {
            int from_row = row - row_step;
            if (from_row + from_col < window_start) break;
            step_score += score_for_moving - penalty(chiron_predictions[from_row], chiron_predictions[from_col]);//  signal[from_row], signal[from_col]);
            double proposed_score = score[from_row][from_col] + step_score;
            if (proposed_score > current_cell) {
              current_cell = proposed_score;
              if (update_come_from) {
                come_from[row][column] = pair<int, int>(from_row, from_col);
              }
            }
          }
        }
        if (row - 1 >= 0 && row_ends[row-1] > column) {
          double step_score = score_for_moving;
          int from_row = row - 1;
          for (int col_step = 1; col_step <= max_steps_back[column] && column - col_step >= row_starts[row]; col_step++) {
            int from_col = column - col_step;
            if (from_row + from_col < window_start) break;
            step_score += score_for_moving - penalty(chiron_predictions[from_row], chiron_predictions[from_col]); //signal[from_row], signal[from_col]);
            double proposed_score = score[from_row][from_col] + step_score;
            if (proposed_score > current_cell) {
              current_cell = proposed_score;
              if (update_come_from) {
                come_from[row][column] = pair<int, int>(from_row, from_col);
              }
            }
          }
        }
      }
    }
  }
  
  fprintf(stderr, "\ndtw done, logging\n");
  
  if (log_filename != "") {
    ofstream out(log_filename.c_str(), ios_base::binary);
    out << signal.size() + 1 << "\n";
    for (int row = 0; row <= static_cast<int>(signal.size()); row++) {
      out << row_starts[row] << " " << row_ends[row] << "\n";
    }
    for (double sig : signal) {
      out << sig << " ";
    }
    out << "\n";
    for (int y = 0; y <= static_cast<int>(signal.size()); y++) {
      for (int x = row_starts[y]; x < row_ends[y]; x++) {
        if (come_from[y][x] == pair<int, int>(-1, -1)) {
          out.put(static_cast<int8_t>(-1));
          out.put(static_cast<int8_t>(-1));
        }
        else {
          if (x - come_from[y][x].second > 127 || y - come_from[y][x].first > 127) {
            fprintf(stderr, "warning: a step too long for encoding into one byte, truncating\n");
          }
          int8_t x_step = min(x - come_from[y][x].second, 127);
          int8_t y_step = min(y - come_from[y][x].first, 127);
          out.put(x_step);
          out.put(y_step);
        }
      }
    }
  }
  
  fprintf(stderr, "getting greedy paths\n");
  vector<vector<pair<int, int>>> result;
  queue<pair<int, int>> unprocessed;
  unprocessed.emplace(min_index_sum, max_index_sum);
  while (!unprocessed.empty()) {
    pair<int, int> interval = unprocessed.front();
    unprocessed.pop();
    
    pair<double, pair<int, int>> best = {0, {-1, -1}};
    for (int index_sum = interval.first; index_sum <= interval.second; index_sum++) {
      if (index_sum % 1000 == 0) fprintf(stderr, "index_sum: %d/%d/%d                    \r", interval.first, index_sum, interval.second);
      for (int row = sum_appears[index_sum]; row <= static_cast<int>(signal.size()) && row_starts[row] + row <= index_sum; row++) {
        int column = index_sum - row;
        score[row][column] = 0;
        pair<int, int> from = come_from[row][column];
        if (from != pair<int,int>(-1, -1)) {
          int from_row = from.first;
          int from_col = from.second;
          double step_score = score_for_moving;
          for (int r = from_row; r < row; r++) {
            for (int c = from_col; c < column; c++) {
              step_score += score_for_moving - penalty(chiron_predictions[r], chiron_predictions[c]); //signal[r], signal[c]);
            }
          }
          score[row][column] = score[from] + step_score;
        }
        best = max(best, {score[row][column], {row, column}});
      }
    }
    if (best.first <= 0) continue;
    
    pair<int, int> position = best.second;
    vector<pair<int, int>> path(1, position);
    int high_sum = position.first + position.second;
    while (come_from[position] != pair<int, int>(-1, -1)) {
      position = come_from[position];
      path.push_back(position);
    }
    int low_sum = position.first + position.second;
    
    for (int index_sum = low_sum+1; index_sum <= high_sum; index_sum++) {
      for (int row = sum_appears[index_sum]; row <= static_cast<int>(signal.size()) && row_starts[row] + row <= index_sum; row++) {
        int column = index_sum - row;
        come_from[row][column] = {-1, -1};
        score[row][column] = 0;
      }
    }
    
    if (interval.first < low_sum) {
      unprocessed.emplace(interval.first, low_sum);
    }
    if (high_sum < interval.second) {
      unprocessed.emplace(high_sum, interval.second);
    }
    //unprocessed.push(interval);
    
    
    reverse(path.begin(), path.end());
    fprintf(stderr, "repeat at: %d, %d - %d, %d\n", path[0].first, path[0].second, path.back().first, path.back().second);
    fprintf(stderr, "removing sum interval %d - %dfrom play\n", low_sum, high_sum);
    if (path.back().first >= path[0].second) {
      result.push_back(path);
    }
    /*else {
      break;
    }*/
  }
  return result;
}
