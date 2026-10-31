// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Which models take a temperature. Current Claude models (Opus 5.5, Sonnet
// 5.5, …) take only their own sampling: a request that carries `temperature`
// is refused with HTTP 400 ("`temperature` is deprecated for this model"). So
// the Workbench's own requests send them none: the dock (useCopilot.js
// copilotModel, through the model port's `temperature: null`, the warmer
// retry included) and the watchdog (useWatchdog.js); the settings drawer
// greys the field out. Every other model keeps the setting (0.2 by default)
// and the warmer retry at ANTI_REPEAT_TEMPERATURE or above (provider/model.js).
// The connection test only lists the models (GET /models): it sends no
// temperature to any model.

import { isAnthropicUrl } from './presets';

const CLAUDE_MODEL = /^claude-/i;

/** True for a Claude model id (`claude-…`) or Anthropic's API (api.anthropic.com): send no temperature. */
export const takesNoTemperature = ({ baseUrl, model } = {}) => (
    (typeof model === 'string' && CLAUDE_MODEL.test(model.trim())) || isAnthropicUrl(baseUrl)
);

/**
 * The temperature a request from these settings ({ baseUrl, model }) asks
 * with: null (none sent) for a model that takes none, else `temperature`
 * (the settings' own by default).
 */
export const temperatureFor = (settings, temperature = settings?.temperature) => (
    takesNoTemperature(settings) ? null : temperature
);
