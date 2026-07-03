import React, { useEffect, useRef, useState } from 'react';
import Markdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import rehypeRaw from 'rehype-raw';
import rehypeSanitize, { defaultSchema } from 'rehype-sanitize';
import { buildChatMessages } from '../prompts/templates';
import { useStreamedReply } from '../useStreamedReply';

// Ask questions about the loaded run; the run context is injected into every
// request, so answers can quote the actual Rwp / convergence numbers. Assistant
// replies render Markdown (GitHub-flavored: tables, lists, code) via
// react-markdown, which builds a React tree and disallows raw HTML — so there is
// still no injection surface; user messages stay plain text. Reasoning models
// stream their chain-of-thought first, shown as a collapsible "Thinking" panel.

// Models often use inline HTML (e.g. <br> for line breaks inside table cells),
// so raw HTML is parsed (rehypeRaw) and then sanitized (rehypeSanitize): safe
// tags like <br>/<sub> render, while scripts, event handlers, and images (no
// external loads from model output) are stripped — the no-injection guarantee
// holds. Order matters: raw → sanitize.
const REMARK_PLUGINS = [remarkGfm];
const REHYPE_PLUGINS = [
    rehypeRaw,
    [rehypeSanitize, { ...defaultSchema, tagNames: (defaultSchema.tagNames || []).filter((tag) => tag !== 'img') }]
];

// Open links in a new tab and wrap tables so wide ones scroll instead of
// breaking the bubble. `node` is dropped so it is not spread onto the DOM.
const MARKDOWN_COMPONENTS = {
    a: ({ node, ...props }) => <a {...props} target="_blank" rel="noreferrer" />,
    table: ({ node, ...props }) => <div className="llm-md-table"><table {...props} /></div>
};

const SUGGESTIONS = [
    'Is there any sign of local distortion?',
    'Which atomic sites are most displaced?',
    'Is the refinement converged?'
];

const ArrowUpIcon = () => (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.2" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
        <path d="M12 19V5M6 11l6-6 6 6" />
    </svg>
);

const StopIcon = () => (
    <svg viewBox="0 0 24 24" fill="currentColor" aria-hidden="true">
        <rect x="7" y="7" width="10" height="10" rx="2" />
    </svg>
);

const UserGlyph = () => (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.9" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
        <circle cx="12" cy="8" r="3.5" />
        <path d="M5 20c0-3.6 3.1-5.5 7-5.5s7 1.9 7 5.5" />
    </svg>
);

// The app's pair-distribution-wave mark, shared with the header and favicon.
const WaveMark = () => (
    <svg viewBox="0 0 100 100" fill="none" aria-hidden="true">
        <path
            d="M14 62 C22 62 27 24 37 24 C47 24 45 66 54 66 C62 66 61 44 69 44 C76 44 78 57 86 57"
            transform="translate(0,5)"
            stroke="currentColor"
            strokeWidth="9"
            strokeLinecap="round"
            strokeLinejoin="round"
        />
    </svg>
);

const ChevronIcon = () => (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
        <path d="M6 9l6 6 6-6" />
    </svg>
);

const formatThinkTime = (ms) => {
    if (!Number.isFinite(ms)) return null;
    const seconds = ms / 1000;
    return seconds < 1 ? '<1s' : `${Math.round(seconds)}s`;
};

// The animated dots shown while the model has started but produced nothing yet
// (a non-reasoning model slow to its first token).
const ThinkingDots = () => (
    <div className="llm-thinking" role="status" aria-label="Thinking">
        <span className="llm-thinking-dot" />
        <span className="llm-thinking-dot" />
        <span className="llm-thinking-dot" />
    </div>
);

// A collapsible chain-of-thought panel: auto-expanded and shimmering while the
// model is still thinking, then collapsed to "Thought for Ns" (re-expandable).
const ReasoningPanel = ({ reasoning, active, reasoningMs }) => {
    const [manualOpen, setManualOpen] = useState(null);
    const open = manualOpen ?? active;
    const bodyRef = useRef(null);

    useEffect(() => {
        if (active && open && bodyRef.current) bodyRef.current.scrollTop = bodyRef.current.scrollHeight;
    }, [reasoning, active, open]);

    const label = active
        ? 'Thinking'
        : (reasoningMs != null ? `Thought for ${formatThinkTime(reasoningMs)}` : 'Reasoning');

    return (
        <div className={`llm-reasoning${open ? ' is-open' : ''}`}>
            <button
                type="button"
                className="llm-reasoning-toggle"
                onClick={() => setManualOpen(!open)}
                aria-expanded={open}
            >
                <span className={`llm-reasoning-label${active ? ' is-active' : ''}`}>{label}</span>
                <span className="llm-reasoning-chevron" aria-hidden="true"><ChevronIcon /></span>
            </button>
            {open && <div className="llm-reasoning-body" ref={bodyRef}>{reasoning}</div>}
        </div>
    );
};

const Message = ({ role, name, content, reasoning, reasoningMs, streaming }) => {
    const isUser = role === 'user';
    const thinkingOnly = streaming && !content;
    return (
        <div className={`llm-message is-${role}`}>
            <span className="llm-message-avatar" aria-hidden="true">
                {isUser ? <UserGlyph /> : <WaveMark />}
            </span>
            <div className="llm-message-body">
                <span className="llm-message-name">{name}</span>
                {reasoning ? (
                    <ReasoningPanel reasoning={reasoning} active={thinkingOnly} reasoningMs={reasoningMs} />
                ) : null}
                {thinkingOnly && !reasoning && <ThinkingDots />}
                {content ? (
                    <div className={`llm-bubble${isUser ? '' : ' is-markdown'}${streaming ? ' is-streaming' : ''}`}>
                        {isUser
                            ? content
                            : (
                                <Markdown
                                    remarkPlugins={REMARK_PLUGINS}
                                    rehypePlugins={REHYPE_PLUGINS}
                                    components={MARKDOWN_COMPONENTS}
                                >
                                    {content}
                                </Markdown>
                            )}
                        {streaming && <span className="llm-caret" aria-hidden="true" />}
                    </div>
                ) : null}
            </div>
        </div>
    );
};

const ChatView = ({ context, settings, connected }) => {
    const [turns, setTurns] = useState([]);
    const [draft, setDraft] = useState('');
    const { text: pending, reasoning: pendingReasoning, reasoningMs: pendingReasoningMs, streaming, error, start, stop } = useStreamedReply();
    const logRef = useRef(null);
    const textareaRef = useRef(null);

    const modelName = settings.model || 'Model';
    const hasConversation = turns.length > 0 || streaming;

    // Keep the newest message in view as turns arrive and tokens (answer or
    // reasoning) stream in.
    useEffect(() => {
        const el = logRef.current;
        if (el) el.scrollTop = el.scrollHeight;
    }, [turns, pending, pendingReasoning]);

    const resizeTextarea = () => {
        const ta = textareaRef.current;
        if (!ta) return;
        ta.style.height = 'auto';
        ta.style.height = `${Math.min(ta.scrollHeight, 160)}px`;
    };

    const send = async () => {
        const question = draft.trim();
        if (!question || streaming || !connected) return;
        setDraft('');
        requestAnimationFrame(resizeTextarea);
        const history = turns;
        setTurns((current) => [...current, { role: 'user', content: question }]);
        const reply = await start({
            baseUrl: settings.baseUrl,
            model: settings.model,
            temperature: settings.temperature,
            apiKey: settings.apiKey,
            messages: buildChatMessages(context, history, question)
        });
        if (reply?.text) {
            setTurns((current) => [...current, {
                role: 'assistant',
                content: reply.text,
                reasoning: reply.reasoning || undefined,
                reasoningMs: reply.reasoningMs
            }]);
        }
    };

    const handleSubmit = (event) => {
        event.preventDefault();
        send();
    };

    const handleKeyDown = (event) => {
        if (event.key === 'Enter' && !event.shiftKey) {
            event.preventDefault();
            send();
        }
    };

    const applySuggestion = (text) => {
        setDraft(text);
        const ta = textareaRef.current;
        if (ta) {
            ta.focus();
            requestAnimationFrame(resizeTextarea);
        }
    };

    return (
        <div className="llm-chat">
            {hasConversation ? (
                <div className="llm-chat-log" ref={logRef}>
                    {turns.map((turn, index) => (
                        <Message
                            key={index}
                            role={turn.role}
                            name={turn.role === 'user' ? 'You' : modelName}
                            content={turn.content}
                            reasoning={turn.reasoning}
                            reasoningMs={turn.reasoningMs}
                        />
                    ))}
                    {streaming && (
                        <Message
                            role="assistant"
                            name={modelName}
                            content={pending}
                            reasoning={pendingReasoning}
                            reasoningMs={pendingReasoningMs}
                            streaming
                        />
                    )}
                </div>
            ) : (
                <div className="llm-chat-empty">
                    <span className="llm-chat-empty-icon" aria-hidden="true"><WaveMark /></span>
                    <p className="llm-chat-empty-title">Ask about this run</p>
                    <p className="llm-chat-empty-hint">
                        {connected
                            ? 'The run’s metrics, symmetry, and convergence history travel with every message, so answers can quote the actual numbers.'
                            : 'Connect a local model — pick one from the selector above — to start a conversation.'}
                    </p>
                    {connected && (
                        <div className="llm-suggestions">
                            {SUGGESTIONS.map((text) => (
                                <button key={text} type="button" className="llm-suggestion" onClick={() => applySuggestion(text)}>
                                    {text}
                                </button>
                            ))}
                        </div>
                    )}
                </div>
            )}

            {error && <div className="llm-connection-status is-error" role="alert">{error}</div>}

            <form className="llm-composer" onSubmit={handleSubmit}>
                <textarea
                    ref={textareaRef}
                    className="llm-composer-input"
                    rows={1}
                    value={draft}
                    placeholder={connected ? 'Message the model about this run…' : 'Connect a model to chat…'}
                    onChange={(event) => {
                        setDraft(event.target.value);
                        resizeTextarea();
                    }}
                    onKeyDown={handleKeyDown}
                    disabled={!connected}
                />
                {streaming ? (
                    <button type="button" className="llm-send is-stop" onClick={stop} aria-label="Stop generating">
                        <StopIcon />
                    </button>
                ) : (
                    <button type="submit" className="llm-send" disabled={!connected || !draft.trim()} aria-label="Send message">
                        <ArrowUpIcon />
                    </button>
                )}
            </form>
        </div>
    );
};

export default ChatView;
