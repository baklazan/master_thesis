#ifndef LOCAL_ALIGNMENT_H
#define LOCAL_ALIGNMENT_H
#include <vector>
#include <string>

std::vector<std::vector<std::pair<int, int>>> 
local_alignment(std::vector<double> signal,
                std::vector<std::vector<float>> chiron_predictions,
                int min_events_distance,
                int max_events_distance,
                double score_for_moving,
                int max_speed_ratio,
                double event_threshold,
                int min_lookahead,
                std::string log_filename="");

#endif
