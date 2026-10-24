// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { buildWatchdogMessages, parseWatchdogReply } from '../prompts/templates';

describe('watchdog prompt round trip', () => {
    it('mentions the stats and the heuristic label', () => {
        const stats = { n_steps: 500, last: 0.87 };
        const messages = buildWatchdogMessages(stats, 'improving', 'stalled');
        expect(messages[1].content).toContain('"n_steps":500');
        expect(messages[1].content).toContain('"improving"');
        expect(messages[1].content).toContain('"stalled"');
    });

    it('names the chi^2 term the history comes from', () => {
        const messages = buildWatchdogMessages({ n_steps: 5 }, 'improving', null, 'X_ray_(R)1');
        expect(messages[1].content).toContain('fit term X_ray_(R)1');
        expect(messages[1].content).toContain('not the total');
        expect(buildWatchdogMessages({ n_steps: 5 }, 'improving', null)[1].content).not.toContain('fit term');
    });

    it('says the stats are the total of every χ² term when they are', () => {
        const messages = buildWatchdogMessages({ n_steps: 5, of: 'total', still_falling: ['F(Q)_1'] }, 'improving', null, 'X_ray_(R)1');
        expect(messages[1].content).toContain('ln of the total chi^2, the sum of every chi^2 term of the .log');
        expect(messages[1].content).not.toContain('fit term X_ray_(R)1');
    });

    it('parses well-formed replies', () => {
        expect(parseWatchdogReply('STATUS: improving — chi dropped to 0.87')).toEqual({
            status: 'improving',
            note: 'chi dropped to 0.87'
        });
        expect(parseWatchdogReply('status: DIVERGING - rising since step 400').status).toBe('diverging');
        expect(parseWatchdogReply('STATUS: converged')).toEqual({ status: 'converged', note: null });
    });

    it('returns null on replies that ignore the format', () => {
        expect(parseWatchdogReply('The modeling run looks fine to me.')).toBeNull();
        expect(parseWatchdogReply('')).toBeNull();
        expect(parseWatchdogReply(null)).toBeNull();
    });
});
