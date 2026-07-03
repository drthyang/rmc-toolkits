import React, { useEffect, useRef, useState } from 'react';
import { buildChatMessages } from '../prompts/templates';
import { useStreamedReply } from '../useStreamedReply';

// Ask questions about the loaded run; the run context is injected into every
// request, so answers can quote the actual Rwp / convergence numbers. Replies
// render as plain pre-wrapped text — no HTML injection surface.

const SUGGESTIONS = [
    'Which dataset fits worst?',
    'Is the refinement converged?',
    'Summarize the fit quality.'
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

const SparkGlyph = () => (
    <svg viewBox="0 0 24 24" fill="currentColor" aria-hidden="true">
        <path d="M12 2l1.7 5.1a5 5 0 0 0 3.2 3.2L22 12l-5.1 1.7a5 5 0 0 0-3.2 3.2L12 22l-1.7-5.1a5 5 0 0 0-3.2-3.2L2 12l5.1-1.7a5 5 0 0 0 3.2-3.2z" />
    </svg>
);

const Message = ({ role, name, content, streaming }) => {
    const isUser = role === 'user';
    return (
        <div className={`llm-message is-${role}`}>
            <span className="llm-message-avatar" aria-hidden="true">
                {isUser ? <UserGlyph /> : <SparkGlyph />}
            </span>
            <div className="llm-message-body">
                <span className="llm-message-name">{name}</span>
                <div className={`llm-bubble${streaming ? ' is-streaming' : ''}`}>
                    {content}
                    {streaming && <span className="llm-caret" aria-hidden="true" />}
                </div>
            </div>
        </div>
    );
};

const ChatView = ({ context, settings, connected }) => {
    const [turns, setTurns] = useState([]);
    const [draft, setDraft] = useState('');
    const { text: pending, streaming, error, start, stop } = useStreamedReply();
    const logRef = useRef(null);
    const textareaRef = useRef(null);

    const modelName = settings.model || 'Model';
    const hasConversation = turns.length > 0 || streaming;

    // Keep the newest message in view as turns arrive and tokens stream in.
    useEffect(() => {
        const el = logRef.current;
        if (el) el.scrollTop = el.scrollHeight;
    }, [turns, pending]);

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
        if (reply) {
            setTurns((current) => [...current, { role: 'assistant', content: reply }]);
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
                        />
                    ))}
                    {streaming && (
                        <Message role="assistant" name={modelName} content={pending} streaming />
                    )}
                </div>
            ) : (
                <div className="llm-chat-empty">
                    <span className="llm-chat-empty-icon" aria-hidden="true"><SparkGlyph /></span>
                    <p className="llm-chat-empty-title">Ask about this run</p>
                    <p className="llm-chat-empty-hint">
                        {connected
                            ? 'The run’s metrics and convergence history travel with every message, so answers can quote the actual numbers.'
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
