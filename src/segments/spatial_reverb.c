#include "../internal.h"
#define PROBE_COUNT 32

struct spatial_reverb_direction{
  float last;
  float *delay;
  uint32_t delay_idx;
  uint32_t delay_length;
  struct biquad_data lpf;
  struct biquad_data apf;
  float gain;
};

struct spatial_reverb_segment_data{
  struct mixed_buffer *in[2];
  struct mixed_buffer *out[2];
  struct spatial_reverb_direction directions[4];
  uint32_t samplerate;
  uint32_t delay_capacity;
  float distance_delay_factor;
  float max_distance_cutoff;

  float probe_angles[PROBE_COUNT];
  float probe_distances[PROBE_COUNT];
  float probe_materials[PROBE_COUNT];
  uint8_t probe_index;
};

static void update_parameters(struct spatial_reverb_segment_data *data, float *distances, float *hit_ratios, float *absorption_rates){
  float distance_delay_factor = data->distance_delay_factor;
  uint32_t delay_capacity = data->delay_capacity;
  uint32_t samplerate = data->samplerate;

  for(uint8_t d=0; d<4; ++d){
    struct spatial_reverb_direction *dir = &data->directions[d];
    uint32_t delay_length = (distance_delay_factor * distances[d]) * samplerate;
    dir->delay_length = MAX(1, MIN(delay_length, delay_capacity));
    dir->gain = hit_ratios[d];

    biquad_lowpass(samplerate, absorption_rates[d]*samplerate, 0.0, &dir->lpf);
    biquad_allpass(samplerate, absorption_rates[d]*samplerate, 1.0, &dir->apf);
  }
}

static inline float gauss(float x){
  const float dev = 0.2;
  return expf(-0.5 * powf(x/dev, 2)) / (dev*sqrt(2*M_PI));
}

static void recompute_parameter(struct spatial_reverb_segment_data *data, float angle, float *distance, float *hit_ratio, float *absorption_rate){
  float distance_sum = 0.0;
  float hit_ratio_sum = 0.0;
  float absorption_sum = 0.0;
  float weight_sum = 0.0;
  float cutoff = data->max_distance_cutoff;

  float *angles = data->probe_angles;
  float *distances = data->probe_distances;
  float *materials = data->probe_materials;
  for(uint8_t i=0; i<PROBE_COUNT; ++i){
    float weight = gauss(angle - angles[i]);
    float distance = distances[i];
    if(distance < cutoff){
      distance_sum += distances[i] * weight;
      hit_ratio_sum += 1 * weight;
      absorption_sum += materials[i] * weight;
    }
    weight_sum += weight;
  }

  if(0.01 < weight_sum){
    float invsum = 1.0 / weight_sum;
    *distance = distance_sum * invsum;
    *hit_ratio = hit_ratio_sum * invsum;
    *absorption_rate = absorption_sum * invsum;
  }else{
    *distance = 0.0;
    *hit_ratio = 0.0;
    *absorption_rate = 0.0;
  }
}

static void recompute_parameters(struct spatial_reverb_segment_data *data){
  float distances[4];
  float hit_ratios[4];
  float absorption_rates[4];

  recompute_parameter(data, M_PI*0.75, &distances[0], &hit_ratios[0], &absorption_rates[0]);
  recompute_parameter(data, M_PI*1.25, &distances[1], &hit_ratios[1], &absorption_rates[1]);
  recompute_parameter(data, M_PI*1.75, &distances[2], &hit_ratios[2], &absorption_rates[2]);
  recompute_parameter(data, M_PI*0.25, &distances[3], &hit_ratios[3], &absorption_rates[3]);
  
  update_parameters(data, distances, hit_ratios, absorption_rates);
}

int spatial_reverb_segment_free(struct mixed_segment *segment){
  struct spatial_reverb_segment_data *data = (struct spatial_reverb_segment_data *)segment->data;
  for(int d=0; d<4; ++d){
    struct spatial_reverb_direction *dir = &data->directions[d];
    if(dir->delay)
      free(dir->delay);
    dir->delay = 0;
  }
  return 1;
}

int spatial_reverb_segment_start(struct mixed_segment *segment){
  struct spatial_reverb_segment_data *data = (struct spatial_reverb_segment_data *)segment->data;
  for(int d=0; d<4; ++d){
    struct spatial_reverb_direction *dir = &data->directions[d];
    dir->last = 0;
    memset(dir->delay, 0, sizeof(float)*data->delay_capacity);
    dir->delay_idx = 0;
    biquad_reset(&dir->lpf);
    biquad_reset(&dir->apf);
  }
  return 1;
}

VECTORIZE int spatial_reverb_segment_mix(struct mixed_segment *segment){
  struct spatial_reverb_segment_data *data = (struct spatial_reverb_segment_data *)segment->data;

  uint32_t samples = UINT32_MAX;
  float *L_in, *R_in, *L_out, *R_out;
  mixed_buffer_request_read(&L_in, &samples, data->in[0]);
  mixed_buffer_request_read(&R_in, &samples, data->in[1]);
  mixed_buffer_request_write(&L_out, &samples, data->out[0]);
  mixed_buffer_request_write(&R_out, &samples, data->out[1]);

  for(uint32_t i=0; i<samples; ++i){
    float L = L_in[i];
    float R = R_in[i];
    // Upmix
    float d_in[4] = {L,L,R,R};
    float d_out[4];
    
    // Mix per direction
    for(uint8_t d=0; d<4; ++d){
      struct spatial_reverb_direction *dir = &data->directions[d];
      float in = d_in[d];
      float *delay = dir->delay;
      uint32_t delay_idx = dir->delay_idx;
      uint32_t delay_length = dir->delay_length;
      
      float sample = dir->last + in;
      float delayed = delay[delay_idx];
      float gained = delayed * dir->gain;
      float lpfd = biquad_sample(gained, &dir->lpf);
      float apfd = biquad_sample(lpfd, &dir->apf);
      float out = apfd;

      delay[delay_idx] = sample;
      dir->delay_idx = (delay_idx+1) % delay_length;
      dir->last = out;
      d_out[d] = out;
    }
    
    // Downmix
    L_out[i] = (d_out[0] + d_out[1])/2.0;
    R_out[i] = (d_out[2] + d_out[3])/2.0;
  }

  mixed_buffer_finish_read(samples, data->in[0]);
  mixed_buffer_finish_read(samples, data->in[1]);
  mixed_buffer_finish_write(samples, data->out[0]);
  mixed_buffer_finish_write(samples, data->out[1]);
  
  return 1;
}

int spatial_reverb_segment_mix_bypass(struct mixed_segment *segment){
  struct spatial_reverb_segment_data *data = (struct spatial_reverb_segment_data *)segment->data;

  mixed_buffer_transfer(data->in[MIXED_LEFT], data->out[MIXED_LEFT]);
  mixed_buffer_transfer(data->in[MIXED_RIGHT], data->out[MIXED_RIGHT]);
  return 1;
}

int spatial_reverb_segment_set_in(uint32_t field, uint32_t location, void *buffer, struct mixed_segment *segment){
  struct spatial_reverb_segment_data *data = (struct spatial_reverb_segment_data *)segment->data;

  switch(field){
  case MIXED_BUFFER:
    if(location < MIXED_LEFT || MIXED_RIGHT < location){
      mixed_err(MIXED_INVALID_LOCATION);
      return 0;
    }
    data->in[location] = (struct mixed_buffer *)buffer;
    return 1;
  default:
    mixed_err(MIXED_INVALID_FIELD);
    return 0;
  }
}

int spatial_reverb_segment_set_out(uint32_t field, uint32_t location, void *buffer, struct mixed_segment *segment){
  struct spatial_reverb_segment_data *data = (struct spatial_reverb_segment_data *)segment->data;

  switch(field){
  case MIXED_BUFFER:
    if(location < MIXED_LEFT || MIXED_RIGHT < location){
      mixed_err(MIXED_INVALID_LOCATION);
      return 0;
    }
    data->out[location] = (struct mixed_buffer *)buffer;
    return 1;
  default:
    mixed_err(MIXED_INVALID_FIELD);
    return 0;
  }
}

int spatial_reverb_segment_get(uint32_t field, void *value, struct mixed_segment *segment){
  struct spatial_reverb_segment_data *data = (struct spatial_reverb_segment_data *)segment->data;
  switch(field){
  case MIXED_SPATIAL_REVERB_DISTANCE_DELAY: *((float *)value) = data->distance_delay_factor; break;
  case MIXED_SPATIAL_REVERB_MAX_DISTANCE_CUTOFF: *((float *)value) = data->max_distance_cutoff; break;
  case MIXED_BYPASS: *((bool *)value) = (segment->mix == spatial_reverb_segment_mix_bypass); break;
  default: mixed_err(MIXED_INVALID_FIELD); return 0;
  }
  return 1;
}

int spatial_reverb_segment_set(uint32_t field, void *value, struct mixed_segment *segment){
  struct spatial_reverb_segment_data *data = (struct spatial_reverb_segment_data *)segment->data;
  switch(field){
  case MIXED_SPATIAL_REVERB_DISTANCE_DELAY: 
    if(*(float *)value < 0.0){
      mixed_err(MIXED_INVALID_VALUE);
      return 0;
    }
    data->distance_delay_factor = *(float *)value;
    break;
  case MIXED_SPATIAL_REVERB_MAX_DISTANCE_CUTOFF:
    data->max_distance_cutoff = *(float *)value;
    break;
  case MIXED_SPATIAL_REVERB_PARAMETERS: 
    update_parameters(data, (float *)value, ((float *)value)+4, ((float *)value)+8);
    break;
  case MIXED_SPATIAL_REVERB_PROBE: {
    float *values = (float *)value;
    uint8_t index = data->probe_index;
    data->probe_angles[index] = fmod(values[0], 2*M_PI);
    data->probe_distances[index] = values[1];
    data->probe_materials[index] = values[2];
    data->probe_index = (index+1) % PROBE_COUNT;
    recompute_parameters(data);
  }
    break;
  case MIXED_BYPASS:
    if(*(bool *)value){
      segment->mix = spatial_reverb_segment_mix_bypass;
    }else{
      segment->mix = spatial_reverb_segment_mix;
    }
    break;
  default:
    mixed_err(MIXED_INVALID_FIELD);
    return 0;
  }
  return 1;
}

int spatial_reverb_segment_info(struct mixed_segment_info *info, struct mixed_segment *segment){
  IGNORE(segment);
  info->name = "spatial_reverb";
  info->description = "Dynamic reverb based on spatial probing.";
  info->flags = MIXED_INPLACE;
  info->min_inputs = 2;
  info->max_inputs = 2;
  info->outputs = 2;
  
  struct mixed_segment_field_info *field = info->fields;
  set_info_field(field++, MIXED_BUFFER,
                 MIXED_BUFFER_POINTER, 1, MIXED_IN | MIXED_OUT | MIXED_SET,
                 "The buffer for audio data attached to the location.");

  set_info_field(field++, MIXED_SPATIAL_REVERB_DISTANCE_DELAY,
                 MIXED_FLOAT, 1, MIXED_SEGMENT | MIXED_SET | MIXED_GET,
                 "How much delay (in seconds) to use per unit of distance.");

  set_info_field(field++, MIXED_SPATIAL_REVERB_MAX_DISTANCE_CUTOFF,
                 MIXED_FLOAT, 1, MIXED_SEGMENT | MIXED_SET | MIXED_GET,
                 "The maximum distance of a probe before it is considered to have missed.");

  set_info_field(field++, MIXED_SPATIAL_REVERB_PARAMETERS,
                 MIXED_FLOAT, 12, MIXED_SEGMENT | MIXED_SET,
                 "Set the spatial reverb parameters.");

  set_info_field(field++, MIXED_SPATIAL_REVERB_PROBE,
                 MIXED_FLOAT, 3, MIXED_SEGMENT | MIXED_SET,
                 "Set a new spatial reverb probe.");

  set_info_field(field++, MIXED_BYPASS,
                 MIXED_BOOL, 1, MIXED_SEGMENT | MIXED_SET | MIXED_GET,
                 "Bypass the segment's processing.");
  
  clear_info_field(field++);
  return 1;
}

MIXED_EXPORT int mixed_make_segment_spatial_reverb(uint32_t samplerate, struct mixed_segment *segment){
  struct spatial_reverb_segment_data *data = calloc(1, sizeof(struct spatial_reverb_segment_data));
  if(!data){
    mixed_err(MIXED_OUT_OF_MEMORY);
    return 0;
  }

  segment->data = data;
  data->samplerate = samplerate;
  data->delay_capacity = samplerate;
  data->distance_delay_factor = 0.0001;
  data->max_distance_cutoff = 1000.0;

  for(int d=0; d<4; ++d){
    struct spatial_reverb_direction *dir = &data->directions[d];
    dir->delay = calloc(data->delay_capacity, sizeof(float));
    dir->delay_length = 1;
    dir->gain = 0.0;
    biquad_lowpass(data->samplerate, data->samplerate, 0, &dir->lpf);
    biquad_allpass(data->samplerate, data->samplerate, 1, &dir->apf);
    if(!dir->delay){
      mixed_err(MIXED_OUT_OF_MEMORY);
      goto cleanup;
    }
  }
  
  segment->free = spatial_reverb_segment_free;
  segment->start = spatial_reverb_segment_start;
  segment->mix = spatial_reverb_segment_mix;
  segment->set_in = spatial_reverb_segment_set_in;
  segment->set_out = spatial_reverb_segment_set_out;
  segment->info = spatial_reverb_segment_info;
  segment->get = spatial_reverb_segment_get;
  segment->set = spatial_reverb_segment_set;
  return 1;

 cleanup:
  spatial_reverb_segment_free(segment);
  segment->data = 0;
  return 0;
}

int __make_spatial_reverb(void *args, struct mixed_segment *segment){
  return mixed_make_segment_spatial_reverb(ARG(uint32_t, 0), segment);
}

REGISTER_SEGMENT(spatial_reverb, __make_spatial_reverb, 1, {
    {.description = "samplerate", .type = MIXED_UINT32}})
