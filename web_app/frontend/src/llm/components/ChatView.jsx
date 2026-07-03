import React, { useState } from 'react';
import { buildChatMessages } from '../prompts/templates';
import { useStreamedReply } from '../useStreamedReply';

// Ask questions about the loaded run; the run context is injected into every
// request, so answers can quote the actual Rwp / convergence numbers. Replies
// render as plain pre-wrapped text — no HTML injection surface.

const ChatView = ({ context, settings, connected }) => {
    const [turns, setTurns] = useState([]);
    const [draft, setDraft] = useState('');
    const { text: pending, streaming, error, start, stop } = useStreamedReply();

    const handleSend = async (event) => {
        event.preventDefault();
        const question = draft.trim();
        if (!question || streaming || !connected) return;
        setDraft('');
        const history = turns;
        setTurns((current) => [...current, { role: 'user', content: question }]);
        const reply = await start({
            baseUrl: settings.baseUrl,
            model: settings.model,
            temperature: settings.temperature,
            messages: buildChatMessages(context, history, question)
        });
        if (reply) {
            setTurns((current) => [...current, { role: 'assistant', content: reply }]);
        }
    };

    return (
        <div className="llm-view">
            {turns.length === 0 && !streaming && (
                <p className="llm-view-hint">
                    {connected
                        ? 'Ask about the loaded run — e.g. “which dataset fits worst?” or “is the refinement converged?”'
                        : 'Connect to a local model in Settings to start a conversation.'}
                </p>
            )}
            <div className="llm-chat-log">
                {turns.map((turn, index) => (
                    <div key={index} className={`llm-chat-turn is-${turn.role}`}>
                        <span className="llm-chat-role">{turn.role === 'user' ? 'You' : settings.model || 'Model'}</span>
                        <div className="llm-reply-text">{turn.content}</div>
                    </div>
                ))}
                {streaming && (
                    <div className="llm-chat-turn is-assistant">
                        <span className="llm-chat-role">{settings.model || 'Model'}</span>
                        <div className="llm-reply-text">{pending || '…'}</div>
                    </div>
                )}
            </div>
            {error && <div className="llm-connection-status is-error" role="alert">{error}</div>}
            <form className="llm-chat-input" onSubmit={handleSend}>
                <input
                    type="text"
                    value={draft}
                    placeholder="Ask about this run…"
                    onChange={(event) => setDraft(event.target.value)}
                    disabled={!connected}
                />
                {streaming ? (
                    <button type="button" className="llm-stop" onClick={stop}>Stop</button>
                ) : (
                    <button type="submit" disabled={!connected || !draft.trim()}>Send</button>
                )}
            </form>
        </div>
    );
};

export default ChatView;
