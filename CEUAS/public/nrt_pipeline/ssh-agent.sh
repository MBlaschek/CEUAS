#
# SSH tools
#
# start a ssh agent
# reconnect to existing ssh-agent

ssh_agentstart() {
    # get active ssh-agent, or launch new
    if [ ! -d "$HOME/.ssh" ]; then
        mkdir -p "$HOME/.ssh" && chmod 700 "$HOME/.ssh"
    fi
    SNAME=$(hostname -s)

    SSH_AGENT_PID=$(ps -fC ssh-agent | grep "ssh-agent -a $HOME/.ssh/ssh-agent-socket-$SNAME" | awk '{print $2}')
    if [ -z "$SSH_AGENT_PID" ]; then
        # kill all ssh-agents (make sure people do not launch too many)
        pgrep -u "$USER" ssh-agent | xargs kill -9 2>/dev/null
        # cleanup old sockets
        rm -f "$HOME/.ssh/ssh-agent-socket-$SNAME" "$HOME/.ssh/ssh-agent-socket" &>/dev/null
        echo "[s] Starting a new ssh-agent"
        eval "$(ssh-agent -a "$HOME/.ssh/ssh-agent-socket-$SNAME")"
    else
        export SSH_AGENT_PID
        export SSH_AUTH_SOCK="$HOME/.ssh/ssh-agent-socket-$SNAME"
        echo "[s] Connecting to existing ssh-agent"
    fi
}

ssh_agentreconnect() {
    for agent in /tmp/ssh-*/agent.* /run/user/"$UID"/keyring/ssh "$HOME"/.ssh/ssh-agent-socket*; do
        if [ -S "$agent" ]; then
            export SSH_AUTH_SOCK="$agent"
            if ssh-add -l &>/dev/null; then
                echo "[s] Found working SSH Agent:"
                ssh-add -l
                return
            fi
        fi
    done
    echo "[s] No ssh-agent found."
    echo "[s] Start a new agent with: ssh_agentstart"
    echo "[s] Note: Please do not start ssh-agent yourself."
}
